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
from tqdm import tqdm


# In[3]:


import matplotlib.pyplot as plt


# # Monte Carlo-minimization
# 
# Li, Zhenqin, and Harold A. Scheraga. "Monte Carlo-minimization approach to the multiple-minima problem in protein folding." Proceedings of the National Academy of Sciences 84, no. 19 (1987): 6611-6615.

# In[4]:


modeller = msm.convert('metenkephalin.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff)


# In[5]:


explorer = oe.Explorer(topology, system, platform='CUDA')


# In[6]:


explorer.set_coordinates(positions)


# In[7]:


exploration_campaign = oe.exploration_campaign.MonteCarloMinimization(explorer)


# In[8]:


exploration_campaign.run(2500, tqdm=True)


# In[12]:


exploration_campaign.acceptance.n_accepted


# In[13]:


exploration_campaign.pes.n_minima


# In[14]:


exploration_campaign.pes.global_minimum_potential_energy


# In[15]:


plt.plot(exploration_campaign.global_minimum_potential_energies._value)


# In[ ]:




