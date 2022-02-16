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


# # Explorer

# In[3]:


topology = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Topology')
forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds, nonbondedMethod=app.NoCutoff)
explorer = oe.Explorer(topology, system, platform='CUDA')


# In[4]:


positions = msm.get('alanine_dipeptide.pdb', coordinates=True)[0]
explorer.set_coordinates(positions)


# In[5]:


explorer.get_potential_energy()


# In[6]:


explorer.get_potential_energy_gradient()


# In[7]:


explorer.get_potential_energy_hessian()


# In[8]:


coordinates = explorer.get_coordinates()


# In[9]:


explorer_2 = explorer.replicate()


# In[10]:


explorer_2


# ## Quenching

# In[ ]:


explorer.set_coordinates(positions)
explorer.quench.l_bfgs()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.quench.fire()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.quench.gradient_descent()
explorer.get_potential_energy()


# ## Moves

# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_atoms_shifts()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_atoms_max_shifts()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_atoms_rsmd()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_atoms_max_rsmd()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_dihedral_shifts()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_dihedral_max_shifts()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_dihedral_rmsd()
explorer.get_potential_energy()


# In[ ]:


explorer.set_coordinates(positions)
explorer.move.random_dihedral_max_rmsd()
explorer.get_potential_energy()


# ## Dynamics

# In[ ]:


explorer.set_coordinates(positions)
explorer.md.langevin(500)
explorer.get_potential_energy()


# ## Distance

# In[ ]:


explorer.set_coordinates(coordinates)
explorer.md.langevin(500)


# In[ ]:


explorer.distance.rmsd(positions)


# In[ ]:


explorer.distance.least_rmsd(positions)


# In[ ]:


explorer.set_coordinates(positions)


# In[ ]:


explorer_2 = explorer.replicate()


# In[ ]:


explorer.md.langevin(500)


# In[ ]:


explorer.distance.rmsd(explorer_2)


# In[ ]:


explorer.distance.least_rmsd(explorer_2)


# In[ ]:




