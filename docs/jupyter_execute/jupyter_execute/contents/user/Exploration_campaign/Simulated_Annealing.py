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


# # Simulated Annealing MonteCarlo
# 
# Kirkpatrick, S., Gelatt, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. science, 220(4598), 671-680.
# 
# ## Goal
# 
# - Global optimization

# In[3]:


modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, #constraints=app.HBonds,
                                implicitSolvent=app.OBC2, soluteDielectric=1.0, solventDielectric=78.5)


# ## Exploration Strategy

# In[ ]:





# ## Exploration Campaign

# In[ ]:




