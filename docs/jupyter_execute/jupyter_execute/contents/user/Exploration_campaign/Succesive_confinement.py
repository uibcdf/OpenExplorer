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


# # Successive confinement
# 
# S. V. Krivov, S. F. Chekmarev and M. Karplus, Phys. Rev. Lett. 88, 3, 2002.
# 

# In[3]:


modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, #constraints=app.HBonds,
                                implicitSolvent=app.OBC2, soluteDielectric=1.0, solventDielectric=78.5)


# In[4]:


explorer = oe.Explorer(topology, system, platform='CUDA')
explorer.set_coordinates(positions)


# In[5]:


exploration = oe.exploration_campaign.SuccessiveConfinement(explorer)


# In[6]:


exploration.run(progress_bar=True, verbose=True)


# Podemos calcular la distancia entre configuraciones para hacer un mapa por proximidad: mds o red con threshold.

# In[ ]:




