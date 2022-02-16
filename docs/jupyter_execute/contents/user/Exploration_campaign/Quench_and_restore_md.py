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
import matplotlib.pyplot as plt


# # Quench and restore
# 
# MD at high temperature with quenching every short periods of time: getting the "hidden" or "inherent structures".
# 
# F. H. Stillinger and T. A. Weber, Phys. Rev. A 25, 978, 1982.    
# F. H. Stillinger and T. A. Weber, J. Phys. Chem. 87, 2833, 1983.    
# F. H. Stillinger and T. A. Weber, Science 225, 983, 1984.    

# In[3]:


## Test system

modeller = msm.convert('alanine_tetrapeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds, nonbondedMethod=app.NoCutoff)

explorer = oe.Explorer(topology, system, platform='CUDA')
explorer.set_coordinates(positions)


# Podemos calcular la distancia entre configuraciones para hacer un mapa por proximidad: mds o red con threshold.

# In[4]:


exploration = oe.exploration_campaign.QuenchAndRestore(explorer)


# In[5]:


exploration.run(2500, progress_bar=True)


# In[ ]:


exploration.pes.n_minima


# In[ ]:


exploration.pes.potential_energy_minima


# In[ ]:


plt.plot(np.array(exploration.time._value)*exploration.time.unit, exploration.trajectory_inherent_structures)
plt.show()


# In[ ]:




