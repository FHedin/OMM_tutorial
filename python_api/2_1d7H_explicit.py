from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# cuton cutoff temperature and pressure chosen for reproducing results from articles of A. Caflisch

prmtop = AmberPrmtopFile('input/1d7h_e.prmtop')
inpcrd = AmberInpcrdFile('input/1d7h_e.inpcrd')

system = prmtop.createSystem(nonbondedMethod=CutoffPeriodic,
		switchDistance=1.0*nanometer,nonbondedCutoff=1.2*nanometer,
                constraints=HBonds)
system.addForce(MonteCarloBarostat(1*bar, 310*kelvin))

integrator = LangevinIntegrator(310*kelvin, 2/picosecond, 0.002*picoseconds)

#platform = Platform.getPlatformByName('CPU')
#platformProperties = {'CpuThreads':'4'}

#platform = Platform.getPlatformByName('OpenCL')
#platformProperties = {'OpenCLPrecision':'single'}

platform = Platform.getPlatformByName('CUDA')
platformProperties = {'CudaPrecision': 'mixed'}

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(inpcrd.positions)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
simulation.minimizeEnergy()

simulation.saveState('output/1d7h_e.minim.state')
simulation.saveCheckpoint('output/1d7h_e.minim.chk')

simulation.context.setVelocitiesToTemperature(310*kelvin)

print('Equilibrating NVT...')

eqsteps = 50000
svfreq = 500
simulation.reporters.append(StateDataReporter(stdout, svfreq, step=True, 
    time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
    volume=True, density=True,
    temperature=True, progress=True,remainingTime=True, speed=True, totalSteps=eqsteps, separator='\t'))
simulation.reporters.append(DCDReporter('output/1d7h_e.npt.dcd', svfreq))
simulation.step(eqsteps)

simulation.saveState('output/1d7h_e.eq_npt_100ps.state')
simulation.saveCheckpoint('output/1d7h_e.eq_npt_100ps.chk')

simulation.currentStep = 0
simulation.context.setTime(0.0)

state = simulation.context.getState(getPositions=True, getVelocities=True,
                 getForces=True, getEnergy=True, getParameters=True, enforcePeriodicBox=True)

#serialize the state to xml for load in parRep later
print('Serializing state : it contains coordinates, velocities, etc. for parRep...')
f = open('output/1d7h_e.State.xml','w')
f.write(XmlSerializer.serialize(state))
f.close()

#serialize the system to xml for load in parRep later
print('Serializing system : it contains all definition of forces for parRep...')
f = open('output/1d7h_e.System.xml','w')
f.write(XmlSerializer.serialize(system))
f.close()

#serialize the integrator to xml for load in parRep later
print('Serializing integrator for parRep...')
f = open('output/1d7h_e.Integrator.xml','w')
f.write(XmlSerializer.serialize(integrator))
f.close()


