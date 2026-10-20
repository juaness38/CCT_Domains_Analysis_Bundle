#!/usr/bin/env python3
"""
Publication-Grade Protein-Peptide Complex Stability Simulation
==============================================================
OpenMM + HMR protocol optimized for kinase-peptide complex simulations.

Protocol (Publication Standard):
1. Energy Minimization (steepest descent + L-BFGS convergence)
2. NVT Equilibration (Gradual heating 0→300K, backbone restraints)
3. NPT Equilibration (Pressure coupling, gradual restraint release)
4. Production MD (Unrestrained, HMR 4fs timestep)

References:
- Amber14SB: Maier et al. (2015) J. Chem. Theory Comput.
- TIP3P: Jorgensen et al. (1983) J. Chem. Phys.
- HMR: Hopkins et al. (2015) J. Chem. Theory Comput.
- LangevinMiddle: Zhang et al. (2019) J. Chem. Phys.

Usage:
    python runcomplex_paper.py --pdb complex_solvated.pdb --output_dir results --ns 200

Author: [Your Name]
Date: 2025
"""

import os
import sys
import argparse
import time
import numpy as np
from datetime import datetime

from openmm import app, unit, CustomExternalForce
import openmm as mm
from openmm.app import PDBFile, PDBxFile, ForceField, Simulation, Modeller
from openmm.app import StateDataReporter, DCDReporter, CheckpointReporter

# Optional: PDBFixer for robust preparation
try:
    import pdbfixer
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False
    print("⚠️  PDBFixer not installed. Using basic preparation.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Publication-grade MD simulation for protein-peptide complexes"
    )
    parser.add_argument('--pdb', required=True, help="Input PDB/CIF file (raw or solvated)")
    parser.add_argument('--output_dir', default="stability_output", help="Output directory")
    parser.add_argument('--replica', type=int, default=1, help="Replica ID for ensemble")
    parser.add_argument('--seed', type=int, default=None, help="Random seed (None=random)")
    parser.add_argument('--ns', type=float, default=100.0, help="Production time (ns)")
    parser.add_argument('--platform', default='CUDA', choices=['CUDA', 'OpenCL', 'CPU'], 
                        help="Compute platform")
    parser.add_argument('--gpu_id', default=None, help="CUDA device index (auto-detect from MPS_GPU_ID env)")
    parser.add_argument('--resume', default=None, help="Resume from checkpoint (.chk)")
    parser.add_argument('--equilibrate_only', action='store_true', 
                        help="Run only equilibration (no production)")
    # Preparation options
    parser.add_argument('--prepare', action='store_true',
                        help="Prepare system: remove waters/ions, add hydrogens, solvate")
    parser.add_argument('--padding', type=float, default=1.2,
                        help="Box padding in nm (default: 1.2)")
    parser.add_argument('--ionic_strength', type=float, default=0.15,
                        help="Ionic strength in M (default: 0.15 M NaCl)")
    parser.add_argument('--ph', type=float, default=7.4,
                        help="pH for protonation states (default: 7.4)")
    return parser.parse_args()


def setup_logging(output_dir, replica):
    """Create log file for simulation metadata."""
    log_path = os.path.join(output_dir, f"replica_{replica}_simulation.log")
    return open(log_path, 'w', encoding='utf-8')


def log_message(log_file, message, also_print=True):
    """Write message to log and optionally print."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"[{timestamp}] {message}"
    log_file.write(formatted + "\n")
    log_file.flush()
    if also_print:
        try:
            print(formatted)
        except UnicodeEncodeError:
            encoding = sys.stdout.encoding or "utf-8"
            safe_message = formatted.encode(encoding, errors="replace").decode(encoding, errors="replace")
            print(safe_message)


def create_restraint_force(topology, positions, k, atom_names=['CA', 'C', 'N'], 
                           exclude_residues=['HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'ZN']):
    """
    Create harmonic positional restraints for specified atoms.
    
    Parameters:
    -----------
    k : float
        Force constant in kJ/mol/nm^2
    atom_names : list
        Atom names to restrain (default: backbone heavy atoms)
    exclude_residues : list
        Residue names to exclude (solvent, ions)
    
    Returns:
    --------
    force : CustomExternalForce
    restrained_atoms : list of atom indices
    """
    # k_restr in kJ/mol/nm^2 - register as float (OpenMM requirement)
    force = CustomExternalForce("k_restr*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k_restr", k)  # k is already in kJ/mol/nm^2 as float
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    
    restrained_atoms = []
    for atom in topology.atoms():
        if atom.name in atom_names and atom.residue.name not in exclude_residues:
            pos = positions[atom.index]
            force.addParticle(atom.index, [pos[0], pos[1], pos[2]])
            restrained_atoms.append(atom.index)
    
    return force, restrained_atoms


def validate_energy(simulation, context_name, log_file, max_energy=1e7):
    """Check that energy is reasonable (not NaN or extremely high)."""
    state = simulation.context.getState(getEnergy=True)
    pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    
    if np.isnan(pe) or np.isinf(pe):
        log_message(log_file, f"❌ ERROR: {context_name} - Energy is NaN/Inf!")
        raise ValueError(f"Invalid energy at {context_name}: {pe}")
    
    if pe > max_energy:
        log_message(log_file, f"⚠️ WARNING: {context_name} - Very high energy: {pe:.2f} kJ/mol")
    else:
        log_message(log_file, f"✓ {context_name} - Energy: {pe:.2f} kJ/mol")
    
    return pe


def run_minimization(simulation, log_file, tolerance=10.0, max_iterations=5000):
    """
    Two-stage minimization: coarse then fine.
    
    Stage 1: Coarse (tolerance=100) - resolve bad contacts
    Stage 2: Fine (tolerance=10) - proper convergence
    """
    log_message(log_file, "🔧 Starting energy minimization...")
    
    # Initial energy
    initial_pe = validate_energy(simulation, "Pre-minimization", log_file, max_energy=1e10)
    
    # Stage 1: Coarse minimization
    log_message(log_file, "   Stage 1: Coarse minimization (tolerance=100 kJ/mol/nm)...")
    simulation.minimizeEnergy(tolerance=100*unit.kilojoule_per_mole/unit.nanometer, 
                               maxIterations=max_iterations)
    validate_energy(simulation, "Post-Stage1", log_file)
    
    # Stage 2: Fine minimization  
    log_message(log_file, f"   Stage 2: Fine minimization (tolerance={tolerance} kJ/mol/nm)...")
    simulation.minimizeEnergy(tolerance=tolerance*unit.kilojoule_per_mole/unit.nanometer, 
                               maxIterations=max_iterations)
    final_pe = validate_energy(simulation, "Post-minimization", log_file)
    
    log_message(log_file, f"   ΔE = {final_pe - initial_pe:.2f} kJ/mol")


def run_nvt_heating(simulation, log_file, target_temp=300, duration_ps=500, 
                    heating_steps=10):
    """
    Gradual NVT heating with temperature ramp.
    
    Parameters:
    -----------
    target_temp : float
        Target temperature in Kelvin
    duration_ps : float
        Total heating duration in picoseconds
    heating_steps : int
        Number of temperature increments
    """
    log_message(log_file, f"🔥 NVT Heating: 0 → {target_temp}K over {duration_ps} ps...")
    
    timestep = simulation.integrator.getStepSize()
    total_steps = int((duration_ps * unit.picosecond) / timestep)
    steps_per_increment = total_steps // heating_steps
    temp_increment = target_temp / heating_steps
    
    # Start from low temperature
    simulation.context.setVelocitiesToTemperature(10 * unit.kelvin)
    
    for i in range(heating_steps):
        current_temp = (i + 1) * temp_increment
        # Update thermostat temperature
        simulation.integrator.setTemperature(current_temp * unit.kelvin)
        simulation.step(steps_per_increment)
        
        # Log progress
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        log_message(log_file, f"   T = {current_temp:.0f}K, PE = {pe:.2f} kJ/mol")
    
    # Final equilibration at target temperature
    log_message(log_file, f"   Equilibrating at {target_temp}K...")
    remaining_steps = total_steps - (steps_per_increment * heating_steps)
    if remaining_steps > 0:
        simulation.step(remaining_steps)
    
    validate_energy(simulation, "Post-NVT heating", log_file)


def run_npt_equilibration(simulation, system, log_file, duration_ps=2000, 
                          restraint_force_index=None, 
                          gradual_release=True):
    """
    NPT equilibration with optional gradual restraint release.
    
    Parameters:
    -----------
    duration_ps : float
        Total NPT duration in picoseconds
    restraint_force_index : int
        Index of restraint force in system
    gradual_release : bool
        If True, gradually reduce restraints to zero
    """
    log_message(log_file, f"💧 NPT Equilibration: {duration_ps} ps at 1 bar...")
    
    timestep = simulation.integrator.getStepSize()
    total_steps = int((duration_ps * unit.picosecond) / timestep)
    
    if gradual_release and restraint_force_index is not None:
        # Gradual release in 4 stages: 100% -> 50% -> 25% -> 10% -> 0%
        stages = [1.0, 0.5, 0.25, 0.1, 0.0]
        steps_per_stage = total_steps // len(stages)
        
        force = system.getForce(restraint_force_index)
        initial_k = 500.0  # kJ/mol/nm^2
        
        for i, fraction in enumerate(stages):
            k = initial_k * fraction
            simulation.context.setParameter("k_restr", k)  # k as float in kJ/mol/nm^2
            simulation.step(steps_per_stage)
            
            state = simulation.context.getState(getEnergy=True)
            pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
            log_message(log_file, 
                f"   Restraints: {fraction*100:.0f}% (k={k:.0f}), PE = {pe:.2f} kJ/mol")
    else:
        simulation.step(total_steps)
    
    validate_energy(simulation, "Post-NPT equilibration", log_file)


def save_equilibrated_structure(simulation, topology, output_path, log_file):
    """Save equilibrated structure as PDB."""
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    
    with open(output_path, 'w') as f:
        PDBFile.writeFile(topology, positions, f)
    
    log_message(log_file, f"💾 Saved equilibrated structure: {output_path}")


def prepare_system(input_pdb, forcefield, log_file, padding=1.2, ionic_strength=0.15, ph=7.4):
    """
    Prepare system from raw PDB: remove crystallographic water/ions, add hydrogens, solvate.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file
    forcefield : ForceField
        OpenMM ForceField object
    padding : float
        Box padding in nm
    ionic_strength : float
        Ionic strength in M
    ph : float
        pH for protonation
        
    Returns:
    --------
    modeller : Modeller
        Prepared system with topology and positions
    """
    log_message(log_file, "\n🧪 SYSTEM PREPARATION")
    log_message(log_file, "=" * 60)
    
    # Load structure
    if input_pdb.endswith('.cif'):
        pdb = PDBxFile(input_pdb)
    else:
        pdb = PDBFile(input_pdb)
    
    log_message(log_file, f"   Loaded: {input_pdb}")
    
    # Use PDBFixer if available (more robust)
    if HAS_PDBFIXER:
        log_message(log_file, "   Using PDBFixer for preparation...")
        fixer = pdbfixer.PDBFixer(filename=input_pdb)
        
        # Find and report missing residues (don't add them for crystallographic structures)
        fixer.findMissingResidues()
        if fixer.missingResidues:
            log_message(log_file, f"   ⚠️ Missing residues detected: {len(fixer.missingResidues)}")
            # Clear missing residues to avoid adding unstructured loops
            fixer.missingResidues = {}
        
        # Find and add missing atoms
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        log_message(log_file, "   ✓ Added missing heavy atoms")
        
        # Remove crystallographic water and heterogens (ligands stay if standard residue)
        fixer.removeHeterogens(keepWater=False)
        log_message(log_file, "   ✓ Removed crystallographic waters and heterogens")
        
        # Add hydrogens at specified pH
        fixer.addMissingHydrogens(ph)
        log_message(log_file, f"   ✓ Added hydrogens at pH {ph}")
        
        modeller = Modeller(fixer.topology, fixer.positions)
    else:
        # Basic preparation without PDBFixer
        log_message(log_file, "   Using basic OpenMM preparation...")
        modeller = Modeller(pdb.topology, pdb.positions)
        
        # Remove water
        modeller.deleteWater()
        log_message(log_file, "   ✓ Removed existing water molecules")
        
        # Add hydrogens
        modeller.addHydrogens(forcefield, pH=ph)
        log_message(log_file, f"   ✓ Added hydrogens at pH {ph}")
    
    # Count atoms before solvation
    n_atoms_protein = modeller.topology.getNumAtoms()
    log_message(log_file, f"   Protein atoms: {n_atoms_protein}")
    
    # Solvate
    log_message(log_file, f"   Solvating with {padding} nm padding, {ionic_strength} M NaCl...")
    modeller.addSolvent(
        forcefield,
        model='tip3p',
        padding=padding * unit.nanometer,
        ionicStrength=ionic_strength * unit.molar,
        positiveIon='Na+',
        negativeIon='Cl-'
    )
    
    n_atoms_total = modeller.topology.getNumAtoms()
    n_waters = (n_atoms_total - n_atoms_protein) // 3  # Approximate
    log_message(log_file, f"   ✓ Added ~{n_waters} water molecules")
    log_message(log_file, f"   Total atoms: {n_atoms_total}")
    
    return modeller


def run_production(simulation, log_file, total_ns, output_dir, run_name, 
                   checkpoint_freq_ns=10):
    """
    Production MD run with periodic checkpointing.
    
    Parameters:
    -----------
    total_ns : float
        Total production time in nanoseconds
    checkpoint_freq_ns : float
        Checkpoint frequency in nanoseconds
    """
    log_message(log_file, f"🚀 Starting Production MD: {total_ns} ns...")
    
    timestep = simulation.integrator.getStepSize()
    total_steps = int((total_ns * unit.nanosecond) / timestep)
    checkpoint_steps = int((checkpoint_freq_ns * unit.nanosecond) / timestep)
    
    chk_path = os.path.join(output_dir, f"{run_name}.chk")
    xml_path = os.path.join(output_dir, f"{run_name}_state.xml")
    
    start_time = time.time()
    steps_done = 0
    
    while steps_done < total_steps:
        steps_to_run = min(checkpoint_steps, total_steps - steps_done)
        simulation.step(steps_to_run)
        steps_done += steps_to_run
        
        # Save checkpoint
        simulation.saveCheckpoint(chk_path)
        simulation.saveState(xml_path)
        
        # Progress report
        elapsed = time.time() - start_time
        ns_done = (steps_done * timestep).value_in_unit(unit.nanosecond)
        speed = ns_done / (elapsed / 86400) if elapsed > 0 else 0  # ns/day
        progress = steps_done / total_steps * 100
        
        state = simulation.context.getState(getEnergy=True)
        pe = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
        
        log_message(log_file, 
            f"   Progress: {progress:.1f}% ({ns_done:.2f} ns) | "
            f"Speed: {speed:.1f} ns/day | PE: {pe:.0f} kJ/mol")
    
    total_time = time.time() - start_time
    log_message(log_file, f"✅ Production complete in {total_time/3600:.2f} hours")


def main():
    args = parse_args()
    
    # MPS Support: Auto-detect GPU from environment
    gpu_id = args.gpu_id if args.gpu_id else os.environ.get('MPS_GPU_ID', '0')
    
    # Setup output directory
    os.makedirs(args.output_dir, exist_ok=True)
    run_name = f"replica_{args.replica}"
    
    # Setup logging
    log_file = setup_logging(args.output_dir, args.replica)
    
    log_message(log_file, "=" * 60)
    log_message(log_file, "KINASE-PEPTIDE COMPLEX MD SIMULATION")
    log_message(log_file, "=" * 60)
    log_message(log_file, f"Input: {args.pdb}")
    log_message(log_file, f"Output: {args.output_dir}")
    log_message(log_file, f"Replica: {args.replica}")
    log_message(log_file, f"Production: {args.ns} ns")
    log_message(log_file, f"Platform: {args.platform} (GPU: {gpu_id})")
    log_message(log_file, f"Seed: {args.seed if args.seed else 'random'}")
    log_message(log_file, f"Prepare system: {args.prepare}")
    log_message(log_file, "=" * 60)

    # =========================================================================
    # FORCE FIELD (needed for both preparation and system creation)
    # Amber14SB protein FF + TIP3P water model
    # =========================================================================
    forcefield = ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')

    # =========================================================================
    # PREPARE OR LOAD STRUCTURE
    # =========================================================================
    if args.prepare:
        # Full preparation pipeline
        modeller = prepare_system(
            args.pdb, forcefield, log_file,
            padding=args.padding,
            ionic_strength=args.ionic_strength,
            ph=args.ph
        )
        topology = modeller.topology
        positions = modeller.positions
        
        # Save prepared structure
        prepared_pdb = os.path.join(args.output_dir, f"{run_name}_prepared.pdb")
        with open(prepared_pdb, 'w') as f:
            PDBFile.writeFile(topology, positions, f)
        log_message(log_file, f"💾 Saved prepared structure: {prepared_pdb}")
    else:
        # Load pre-solvated structure
        log_message(log_file, f"📂 Loading structure: {args.pdb}")
        
        if args.pdb.endswith('.cif'):
            pdb = PDBxFile(args.pdb)
        else:
            pdb = PDBFile(args.pdb)
        
        topology = pdb.topology
        positions = pdb.positions
    
    # Count atoms
    n_atoms = sum(1 for _ in topology.atoms())
    n_residues = sum(1 for _ in topology.residues())
    n_chains = sum(1 for _ in topology.chains())
    log_message(log_file, f"   Atoms: {n_atoms}, Residues: {n_residues}, Chains: {n_chains}")

    # =========================================================================
    # CREATE SYSTEM
    # =========================================================================
    log_message(log_file, "\n📦 Creating System...")
    log_message(log_file, "   Force Field: Amber14SB (protein.ff14SB) + TIP3P")
    log_message(log_file, "   HMR: Enabled (4 amu hydrogen mass)")
    
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True,
        hydrogenMass=4.0 * unit.amu,  # HMR for 4fs timestep
        ewaldErrorTolerance=5e-5      # Tight PME tolerance
    )
    
    log_message(log_file, f"   PME cutoff: 1.0 nm")
    log_message(log_file, f"   Constraints: HBonds")
    log_message(log_file, f"   Number of forces: {system.getNumForces()}")

    # Add CMMotionRemover to prevent center-of-mass drift in long simulations
    # This is standard practice for production MD (removes COM motion every 100 steps)
    system.addForce(mm.CMMotionRemover(100))
    log_message(log_file, "   CMMotionRemover: Enabled (every 100 steps)")

    # =========================================================================
    # ADD RESTRAINTS (for equilibration)
    # =========================================================================
    log_message(log_file, "\n🔒 Adding backbone restraints (k=500 kJ/mol/nm²)...")
    restraint_force, restrained_atoms = create_restraint_force(
        topology, positions, k=500.0
    )
    restraint_force_index = system.addForce(restraint_force)
    log_message(log_file, f"   Restrained {len(restrained_atoms)} backbone atoms (CA, C, N)")

    # =========================================================================
    # INTEGRATOR & PLATFORM
    # =========================================================================
    # Use 2fs for equilibration (safer), switch to 4fs for production
    temperature = 300 * unit.kelvin
    friction = 1.0 / unit.picosecond
    timestep_eq = 2.0 * unit.femtoseconds
    
    integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep_eq)
    integrator.setConstraintTolerance(1e-6)  # Tight tolerance for stability
    
    if args.seed:
        integrator.setRandomNumberSeed(args.seed)
        log_message(log_file, f"\n🎲 Random seed: {args.seed}")
    
    # Platform setup (MPS compatible)
    if args.platform == 'CUDA':
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'Precision': 'mixed', 'DeviceIndex': gpu_id}
    elif args.platform == 'OpenCL':
        platform = mm.Platform.getPlatformByName('OpenCL')
        properties = {'Precision': 'mixed'}
    else:
        platform = mm.Platform.getPlatformByName('CPU')
        properties = {}
    
    log_message(log_file, f"\n⚙️ Platform: {platform.getName()}")

    # =========================================================================
    # CREATE SIMULATION
    # =========================================================================
    simulation = Simulation(topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    
    # Check for resume
    if args.resume and os.path.exists(args.resume):
        log_message(log_file, f"📥 Resuming from checkpoint: {args.resume}")
        simulation.loadCheckpoint(args.resume)
    else:
        # =====================================================================
        # PHASE 1: MINIMIZATION
        # =====================================================================
        log_message(log_file, "\n" + "=" * 60)
        log_message(log_file, "PHASE 1: ENERGY MINIMIZATION")
        log_message(log_file, "=" * 60)
        
        run_minimization(simulation, log_file, tolerance=10.0, max_iterations=10000)
        
        # =====================================================================
        # PHASE 2: NVT HEATING
        # =====================================================================
        log_message(log_file, "\n" + "=" * 60)
        log_message(log_file, "PHASE 2: NVT HEATING (0 → 300K)")
        log_message(log_file, "=" * 60)
        
        run_nvt_heating(simulation, log_file, target_temp=300, duration_ps=500, 
                       heating_steps=10)
        
        # =====================================================================
        # PHASE 3: NPT EQUILIBRATION
        # =====================================================================
        log_message(log_file, "\n" + "=" * 60)
        log_message(log_file, "PHASE 3: NPT EQUILIBRATION (1 bar)")
        log_message(log_file, "=" * 60)
        
        # Add barostat NOW (not before NVT)
        pressure = 1.0 * unit.bar
        barostat = mm.MonteCarloBarostat(pressure, temperature, 25)
        barostat_index = system.addForce(barostat)
        simulation.context.reinitialize(preserveState=True)
        
        log_message(log_file, "   Barostat: MonteCarloBarostat (25 step interval)")
        
        run_npt_equilibration(simulation, system, log_file, duration_ps=2000,
                             restraint_force_index=restraint_force_index,
                             gradual_release=True)
        
        # Save equilibrated structure
        eq_pdb_path = os.path.join(args.output_dir, f"{run_name}_equilibrated.pdb")
        save_equilibrated_structure(simulation, topology, eq_pdb_path, log_file)
        
        # Save equilibrated state
        eq_state_path = os.path.join(args.output_dir, f"{run_name}_equilibrated.xml")
        simulation.saveState(eq_state_path)
        log_message(log_file, f"💾 Saved equilibrated state: {eq_state_path}")
    
    if args.equilibrate_only:
        log_message(log_file, "\n✅ Equilibration complete. Exiting (--equilibrate_only)")
        log_file.close()
        return

    # =========================================================================
    # PHASE 4: PRODUCTION
    # =========================================================================
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "PHASE 4: PRODUCTION MD")
    log_message(log_file, "=" * 60)
    
    # Switch to 4fs timestep for HMR
    timestep_prod = 4.0 * unit.femtoseconds
    integrator.setStepSize(timestep_prod)
    log_message(log_file, f"   Timestep: 4 fs (HMR enabled)")
    
    # Setup reporters
    dcd_path = os.path.join(args.output_dir, f"{run_name}.dcd")
    csv_path = os.path.join(args.output_dir, f"{run_name}_energy.csv")
    
    # Save every 10 ps (2500 steps @ 4fs)
    report_freq = 2500
    log_message(log_file, f"   DCD output: {dcd_path} (every 10 ps)")
    log_message(log_file, f"   Energy log: {csv_path}")
    
    simulation.reporters.append(DCDReporter(dcd_path, report_freq))
    simulation.reporters.append(StateDataReporter(
        csv_path, report_freq,
        step=True, time=True, 
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, density=True,
        speed=True
    ))
    simulation.reporters.append(StateDataReporter(
        sys.stdout, report_freq * 10,  # Print every 100 ps
        step=True, time=True, temperature=True, speed=True
    ))
    
    # Run production
    run_production(simulation, log_file, args.ns, args.output_dir, run_name,
                  checkpoint_freq_ns=10)
    
    # Final save
    final_state = os.path.join(args.output_dir, f"{run_name}_final.xml")
    simulation.saveState(final_state)
    
    log_message(log_file, "\n" + "=" * 60)
    log_message(log_file, "SIMULATION COMPLETE")
    log_message(log_file, "=" * 60)
    
    # Summary
    log_message(log_file, "\n📊 Output files:")
    log_message(log_file, f"   - Trajectory: {dcd_path}")
    log_message(log_file, f"   - Energy data: {csv_path}")
    log_message(log_file, f"   - Final state: {final_state}")
    log_message(log_file, f"   - Equilibrated PDB: {eq_pdb_path}")
    
    log_file.close()


if __name__ == "__main__":
    main()
