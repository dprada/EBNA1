import warnings
warnings.filterwarnings("ignore")

import os
import molmodmt as m3t
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from openmmtools.integrators import LangevinIntegrator

# Loading initial PDB file

molmod_modeller = m3t.load('system_init.pdb','openmm.Modeller')

# System 

topology = m3t.convert(molmod_modeller, 'openmm.Topology')
positions = m3t.get(molmod_modeller, coordinates=True)
forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds)

# Thermodynamic state (meaningless here but necessary)

temperature = 0.0*unit.kelvin
pressure    = None

# Integrator

friction   = 1.0/unit.picosecond
step_size  = 2.0*unit.femtoseconds
integrator = LangevinIntegrator(temperature, friction, step_size)
integrator.setConstraintTolerance(0.00001)

# Platform

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

# Simulation

simulation = app.Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(0.0*unit.kelvin)

# Reporting Initial Potential Energy

state = simulation.context.getState(getEnergy=True)
potential_energy_before_minimization = state.getPotentialEnergy()

# Running minimization

simulation.minimizeEnergy()

# Reporting Final Potential Energy

state = simulation.context.getState(getEnergy=True)
potential_energy_after_minimization = state.getPotentialEnergy()

# Output as PDB file

m3t.convert(simulation, 'system_minimized.pdb')

# Print out info
print("Potential Energy before minimization: {}".format(potential_energy_before_minimization))
print("Potential Energy after minimization: {}".format(potential_energy_after_minimization))

