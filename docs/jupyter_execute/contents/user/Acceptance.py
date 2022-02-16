#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


import molsysmt as msm
import openexplorer as oe
import numpy as np
from simtk import unit
from simtk.openmm import app


# # Acceptance
# 
# ## Metropolis-Hastings

# In[3]:


modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds, nonbondedMethod=app.NoCutoff)

explorer = oe.Explorer(topology, system, platform='CUDA')

explorer.set_coordinates(positions)


# In[4]:


msm.get_dihedral_angles(explorer, dihedral_angle='all')


# In[5]:


explorer.get_potential_energy()


# In[6]:


explorer.move.dihedral_shifts()


# In[7]:


msm.get_dihedral_angles(explorer, dihedral_angle='all')


# In[8]:


explorer.get_potential_energy()


# In[9]:


explorer.acceptance.metropolis_hastings(previous_coordinates=positions)


# In[10]:


explorer.acceptance.metropolis_hastings.accepted


# In[11]:


explorer.acceptance.metropolis_hastings.n_tries


# In[12]:


explorer.acceptance.metropolis_hastings.n_accepted


# In[13]:


explorer.acceptance.metropolis_hastings.random


# In[14]:


explorer.acceptance.metropolis_hastings.weight


# In[15]:


explorer.acceptance.metropolis_hastings.potential_energy


# In[16]:


explorer.acceptance.metropolis_hastings.previous_potential_energy


# In[17]:


msm.get_dihedral_angles(explorer, dihedral_angle='all')


# In[18]:


explorer.acceptance.metropolis_hastings._kT


# In[41]:


acceptor=explorer.acceptance.metropolis_hastings
acceptor.run(previous_potential_energy=1.0*unit.kilojoules_per_mole,
             potential_energy=2.5*unit.kilojoules_per_mole,
             update_explorer=False)

print(acceptor.weight, acceptor.random, acceptor.accepted)


# In[ ]:




