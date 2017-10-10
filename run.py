#!/usr/bin/env python2.7

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout

# Input Files
print('Parsing PDB file...')
pdb = PDBFile('ala2.pdb')

print('Parsing FF file(s)...')
forcefield = ForceField('amber99sbildn.xml')

# System Configuration
nonbondedMethod = CutoffNonPeriodic
nonbondedCutoff = 1.4*nanometers
constraints = HBonds
constraintTolerance = 0.00001

# Prepare the Simulation
print('Building system...')
topology = pdb.topology
positions = pdb.positions
system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod, nonbondedCutoff=nonbondedCutoff,
    constraints=constraints)

# Integration Options
dt = 0.002*picoseconds
temperature = 300.00*kelvin
friction = 2/picosecond
integrator = LangevinIntegrator(temperature, friction, dt)
integrator.setConstraintTolerance(constraintTolerance)

# do minimization, perform equilibration, then save the state of the simulation
equilibrationSteps = 25000
platform = Platform.getPlatformByName('CPU')
platformProperties = {'CpuThreads':'1'}

simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

# Minimize and Equilibrate
print('Performing energy minimization...')
simulation.minimizeEnergy()
print('Equilibrating...')
simulation.context.setVelocitiesToTemperature(temperature)
simulation.step(equilibrationSteps)

simulation.currentStep = 0
simulation.context.setTime(0.0)

simulation.saveState('eq.state')
simulation.saveCheckpoint('eq.chk')

simsteps = 50e6 # 100 ns
svfreq=500

simulation.reporters.append(StateDataReporter(stdout, svfreq, step=True, 
    time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    temperature=True, progress=True,remainingTime=True, speed=True, totalSteps=simsteps, separator='\t'))
simulation.reporters.append(DCDReporter('sim_100ns.dcd', svfreq))
simulation.step(simsteps)

simulation.saveState('sim_100ns.state')
simulation.saveCheckpoint('sim_100ns.chk')


