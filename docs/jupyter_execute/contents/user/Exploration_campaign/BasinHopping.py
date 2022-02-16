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


# # BasinHopping

# In[3]:


modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds, nonbondedMethod=app.NoCutoff)

explorer = oe.Explorer(topology, system, platform='CUDA')


# In[4]:


explorer.set_coordinates(positions)


# In[5]:


movement = oe.tools.move.DihedralShifts()


# In[5]:


movement = oe.movements.RandomDihedralMaxShifts(stepsize=30*unit.degrees)


# In[ ]:





# In[ ]:


exploration = oe.exploration_campaign.BasinHopping(explorer, movement, temperature=300.0*unit.kelvin)


# In[ ]:


exploration.run(100)


# In[ ]:


exploration.n_tries


# In[ ]:


exploration.n_acceptances


# In[ ]:


exploration.reset(temperature=900.0*unit.kelvin)


# In[ ]:


exploration.run(50000)


# In[ ]:


exploration.n_tries


# In[ ]:


exploration.n_acceptances


# In[ ]:


np.exp(1./3.0)


# In[ ]:




