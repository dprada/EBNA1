import warnings
warnings.filterwarnings("ignore")

import os
from time import time as realtime
from time import asctime, localtime
import numpy as np
import molmodmt as m3t
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from openmmtools.integrators import LangevinIntegrator


start_realtime = realtime()
print("")
print("Start:",asctime(localtime()))

#### Loading PDB

pdb = m3t.convert('system_minimized.pdb', 'openmm.PDBFile')

#### System

topology = m3t.convert(pdb, 'openmm.Topology')
forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
system = forcefield.createSystem(topology,
                                 nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.2*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)

#### Custom Forces

#### Thermodynamic State

kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
temperature = 300.0*unit.kelvin
pressure    = None

#### Integrator

friction   = 1.0/unit.picosecond
step_size  = 2.0*unit.femtoseconds
integrator = LangevinIntegrator(temperature, friction, step_size)
integrator.setConstraintTolerance(0.00001)

#### Platform

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

#### Simulation

simulation = app.Simulation(topology, system, integrator, platform, properties)

#### Initial Conditions

positions = m3t.get(pdb, coordinates=True)
velocities = None # Set according to Maxwell Boltzmann distribution
box_vectors = None # box_vectors in topology coming from pdb
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(0.0*unit.kelvin)
#simulation.context.setPeriodicBoxVectors(box_vectors) ## r_boxes[0] * nanometer,r_boxes[1] * nanometer, r_boxes[2] * nanometer

#### Iterations Parameters

time_simulation = 100.0 * unit.picoseconds
time_saving = 10.0 * unit.picoseconds
time_verbose = 10.0 * unit.picoseconds
time_checkpoint = 10.0 * unit.picoseconds

number_iterations = int(time_simulation/time_saving)
elapsed_iterations_verbose = int(time_verbose/time_saving)
elapsed_iterations_checkpoint = int(time_checkpoint/time_saving)
steps_per_iteration = int(time_saving/step_size)
total_simulation_steps = number_iterations*steps_per_iteration

#### Reporters

# Observables Stored

net_mass, n_degrees_of_freedom = m3t.get(system, net_mass=True, n_degrees_of_freedom=True)
number_states_saved = number_iterations+1
data = dict()
data['time'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.picoseconds)
data['potential'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.kilocalories_per_mole)
data['kinetic'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.kilocalories_per_mole)
data['volume'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.nanometers**3)
data['density'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.gram / unit.centimeters**3)
data['kinetic_temperature'] = unit.Quantity(np.zeros([number_states_saved], np.float64), unit.kelvin)

# Simulation Status

def printout_status(iteration, number_iterations, time, kinetic_temperature, volume):

    progression = 100.0*(interation/number_iterations)
    print("Progress: {:6.2f}, Iteration: {:d}, Time: {:.2f} ps, Kinetic Temperature: {:.3f} K,\
          Volume: {:.3f} nm^3".format(progression, iteration,
                                      time.value_in_units(unit.picoseconds),
                                     kinetic_temperature.value_in_units(unit.kelvin),
                                     volume.value_in_units(unit.nanometers**3)))

#### CheckPoint and Finnal State

def dump_state(simulation, filename):

    state = simulation.context.getState(getPositions=True, getVelocities=True)
    time = state.getTime() / unit.picoseconds
    positions = state.getPositions() / unit.nanometers
    velocities = state.getVelocities() / unit.nanometers*unit.picoseconds
    box_vectors = state.getPeriodicBoxVectors() / unit.nanometers

    with open(os.path.join(filename), 'wb') as f:
        pickle.dump((time, positions, velocities, box_vectors), f)

def save_checkpoint_state(simulation):
    return dump_state(simulation, "checkpoint.pkl")

def save_finnal_state(simulation):
    return dump_state(simulation, "restart.pkl")

#### Reporting Initial State

state = simulation.context.getState(getEnergy=True)
time = state.getTime()
potential_energy = state.getPotentialEnergy()
kinetic_energy = state.getKineticEnergy()
volume = state.getPeriodicBoxVolume()
density = (net_mass / volume).in_units_of(unit.gram / unit.centimeter**3)
kinetic_temperature = (2.0 * kinetic_energy / kB / n_degrees_of_freedom).in_units_of(unit.kelvin)
data['time'][0] = time
data['potential'][0] = potential_energy
data['kinetic'][0] = kinetic_energy
data['volume'][0] = volume
data['density'][0] = density
data['kinetic_temperature'][0] = kinetic_temperature

#### Running Simulation

start_simulation_realtime = realtime()

for iteration in range(1, number_iterations+1):
    integrator.step(steps_per_iteration)
    state = simulation.context.getState(getEnergy=True)
    time = state.getTime()
    potential_energy = state.getPotentialEnergy()
    kinetic_energy = state.getKineticEnergy()
    volume = state.getPeriodicBoxVolume()
    density = (net_mass / volume).in_units_of(unit.gram / unit.centimeter**3)
    kinetic_temperature = (2.0 * kinetic_energy / kB / n_degrees_of_freedom).in_units_of(unit.kelvin)
    data['time'][iteration]=time
    data['potential'] = potential_energy
    data['kinetic'] = kinetic_energy
    data['volume'] = volume
    data['density'] = density
    data['kinetic_temperature'] = kinetic_temperature
    if (iteration%elapsed_iterations_verbose)==0:
        printout_status(iteration, number_iterations, time, kinetic_temperature, volume)
    if (iteration%elapsed_iterations_checkpoing)==0:
        save_checkpoint_state(simulation)

end_simulation_realtime = realtime()

#### Saving Finnal State

save_finnal_state(simulation)

#### Summary


end_realtime = realtime()
preparation_elapsed_realtime = (start_simulation_realtime - start_realtime)*unit.seconds
simulation_elapsed_realtime = (end_simulation_realtime - start_simulation_realtime)*unit.seconds
total_elapsed_realtime = (end_realtime - start_realtime)*unit.seconds

performance = (time_simulation/unit.nanoseconds) / (simulation_elapsed_realtime/unit.hours)

print("End:",asctime(localtime()))
print("")
print("************************")
print("")
print("Total time: {}".format(total_elapsed_realtime))
print("Preparation time: {}".format(preparation_elapsed_realtime))
print("Simulation time: {}".format(simulation_elapsed_realtime))
print("Performance: {} ns/h", simulation_elapsed_realtime)

