#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import molsysmt as msm
import numpy as np
from simtk import unit
from simtk.openmm import app
import simtk.openmm as mm
from mdtraj.reporters import HDF5Reporter
from sys import stdout
import matplotlib.pyplot as plt


# In[3]:


modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml') # 'amber10_obc.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff,
                                 implicitSolvent=app.OBC2, soluteDielectric=1.0, solventDielectric=78.5)


# In[4]:


temperature = 500*unit.kelvin
integration_timestep = 2.0*unit.femtosecond
collision_rate = 1.0/unit.picosecond

integrator = mm.LangevinIntegrator(temperature, collision_rate, integration_timestep)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed'}
simulation = app.Simulation(topology, system, integrator, platform, properties)
simulation.context.setPositions(positions)


# In[5]:


simulation.minimizeEnergy()


# In[6]:


simulation.context.setVelocitiesToTemperature(temperature)


# In[7]:


md_time = 1.000 * unit.nanoseconds
saving_time_interval = 1.0 * unit.picoseconds
log_time_interval = 50.0 * unit.picoseconds


# In[8]:


md_steps = np.rint(md_time/integration_timestep).astype(int)
saving_steps_interval = np.rint(saving_time_interval/integration_timestep).astype(int)
log_steps_interval = np.rint(log_time_interval/integration_timestep).astype(int)

reporter_log_stdout = app.StateDataReporter(stdout, log_steps_interval, step=True, time=True,
                                            totalEnergy=True, temperature=True,
                                            progress=True, remainingTime=True,
                                            speed=True, totalSteps=md_steps, separator='\t')


reporter_logfile = app.StateDataReporter('traj.log', log_steps_interval, step=True, time=True,
                                            potentialEnergy=True, kineticEnergy=True,
                                            totalEnergy=True, temperature=True,
                                            progress=True, remainingTime=True,
                                            speed=True, totalSteps=md_steps, separator='\t')

reporter_h5file = HDF5Reporter('traj.h5', saving_steps_interval, coordinates=True, time=True,
                               potentialEnergy=True, kineticEnergy=True, temperature=True)

simulation.reporters.append(reporter_log_stdout)
simulation.reporters.append(reporter_logfile)
simulation.reporters.append(reporter_h5file)

simulation.step(md_steps)

reporter_h5file.close()


# In[9]:


molecular_system = msm.convert('traj.h5', to_form='molsysmt.MolSys')


# In[10]:


phi_chains, psi_chains, phi_angles, psi_angles = msm.ramachandran_angles(molecular_system)


# In[11]:


plt.plot(phi_angles)


# In[12]:


plt.plot(psi_angles)


# In[13]:


import seaborn as sns

ax = sns.kdeplot(phi_angles[:,0], psi_angles[:,0], shade=True)
ax.set_xlim(-180.0,180.0)
ax.set_ylim(-180.0,180.0)
ax.set_xticks([-180.0, -90.0, 0.0, 90.0, 180.0])
ax.set_yticks([-180.0, -90.0, 0.0, 90.0, 180.0])
ax.set_xlabel('$\phi_1$')
ax.set_ylabel('$\psi_1$')
ax.set_aspect('equal')


# In[15]:


bb=msm.least_rmsd(molecular_system)


# In[16]:


plt.plot(bb)


# In[ ]:




