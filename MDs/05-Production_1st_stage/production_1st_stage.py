import warnings
warnings.filterwarnings("ignore")

from sys import stdout as _stdout
from sys import exit
from time import time as realtime
from pytools_uibcdf.Time import formatted_elapsed_time, formatted_local_time
import molmodmt as m3t
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
from mdtraj.reporters import HDF5Reporter
from openmmtools.integrators import LangevinIntegrator

#### Log file

logfile = open('logfile.txt','w')
#logfile = _stdout

start_realtime = realtime()
logfile.write("\n")
logfile.write("Start: "+formatted_local_time()+"\n")
logfile.write("\n")

#### Loading PDB

pdb = m3t.convert('system_equilibrated_NVT.pdb', 'openmm.PDBFile')

#### System

topology = m3t.convert(pdb, 'openmm.Topology')
forcefield = app.ForceField('amber99sbildn.xml','tip3p.xml')
system = forcefield.createSystem(topology,
                                 nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.2*unit.nanometers,
                                 constraints=app.HBonds,
                                 rigidWater=True,
                                 ewaldErrorTolerance=0.0005)

#### Thermodynamic State

kB = unit.BOLTZMANN_CONSTANT_kB * unit.AVOGADRO_CONSTANT_NA
temperature = 300.0*unit.kelvin
pressure    = 1.0*unit.atmosphere

#### Integrator

friction   = 1.0/unit.picosecond
step_size  = 2.0*unit.femtoseconds
integrator = LangevinIntegrator(temperature, friction, step_size)
integrator.setConstraintTolerance(0.00001)

#### Barostat

barostat_interval = 25
barostat = mm.MonteCarloBarostat(pressure, temperature, barostat_interval)
system.addForce(barostat)

#### Platform

platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}

#### Simulation

simulation = app.Simulation(topology, system, integrator, platform, properties)

#### Initial Conditions

positions = m3t.get(pdb, coordinates=True)
simulation.context.setPositions(positions)
simulation.context.setVelocitiesToTemperature(temperature)

#### Iterations Parameters

steps_simulation = 10000
steps_interval_saving = 500
steps_interval_verbose = 2500
steps_interval_checkpoint = 2500

time_simulation = (steps_simulation*step_size).in_units_of(unit.picoseconds)
time_saving = (steps_interval_saving*step_size).in_units_of(unit.picoseconds)
time_verbose = (steps_interval_verbose*step_size).in_units_of(unit.picoseconds)
time_checkpoint = (steps_interval_checkpoint*step_size).in_units_of(unit.picoseconds)

logfile.write("\n")
logfile.write("Step size: {}\n".format(step_size))
logfile.write("Simulation time: {} ({} steps)\n".format(time_simulation , steps_simulation))
logfile.write("Saving time: {} ({} steps)\n".format(time_saving , steps_interval_saving))
logfile.write("Verbose time: {} ({} steps)\n".format(time_verbose , steps_interval_verbose))
logfile.write("Checkpoint time: {} ({} steps)\n".format(time_checkpoint , steps_interval_checkpoint))
logfile.write("\n")


#### Reporters

# Logfile

simulation.reporters.append(app.StateDataReporter(logfile, reportInterval=steps_interval_verbose,
                                                  progress=True, speed=True, step=True, time=True,
                                                  potentialEnergy=True, temperature=True,
                                                  volume=True, totalSteps=steps_simulation,
                                                  separator=", "))

# Observables

simulation.reporters.append(HDF5Reporter('traj_1st_stage.h5', reportInterval=steps_interval_saving,
                                         coordinates=True, time=True, cell=True,
                                         potentialEnergy=True, kineticEnergy=True,
                                         temperature=True))

# Checkpoints

simulation.reporters.append(app.CheckpointReporter('checkpnt.chk', steps_interval_checkpoint))

#### Running Simulation

start_simulation_realtime = realtime()

simulation.step(steps_simulation)

end_simulation_realtime = realtime()

#### Saving Finnal State

simulation.saveState('finnal_state.xml')
simulation.saveCheckpoint('finnal_state.chk')
m3t.convert(simulation,'finnal_positions.pdb')

#### Summary

end_realtime = realtime()
preparation_elapsed_realtime = (start_simulation_realtime - start_realtime)*unit.seconds
simulation_elapsed_realtime = (end_simulation_realtime - start_simulation_realtime)*unit.seconds
total_elapsed_realtime = (end_realtime - start_realtime)*unit.seconds

performance = 24 * (steps_simulation*step_size/unit.nanoseconds) / (simulation_elapsed_realtime/unit.hours)

logfile.write("\n")
logfile.write("End: "+formatted_local_time()+"\n")
logfile.write("\n")
logfile.write("****SUMMARY****\n")
logfile.write("\n")
logfile.write("Total time: "+formatted_elapsed_time(total_elapsed_realtime/unit.seconds)+"\n")
logfile.write("Preparation time: "+formatted_elapsed_time(preparation_elapsed_realtime/unit.seconds)+"\n")
logfile.write("Simulation time: "+formatted_elapsed_time(simulation_elapsed_realtime/unit.seconds)+"\n")
logfile.write("\n")
logfile.write("Simulation Performance: {:.3f} ns/day".format(performance)+"\n")
logfile.write("\n")

#### Closing reporters
logfile.close()
simulation.reporters[1].close()

