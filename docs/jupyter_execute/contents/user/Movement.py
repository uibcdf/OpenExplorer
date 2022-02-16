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


# # Movement

# In[3]:


topology = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Topology')
forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')
system = forcefield.createSystem(topology, constraints=app.HBonds, nonbondedMethod=app.NoCutoff)
explorer = oe.Explorer(topology, system, platform='CUDA')


# In[4]:


initial_positions = msm.get('alanine_dipeptide.pdb', coordinates=True)[0]
explorer.set_coordinates(initial_positions)


# In[5]:


explorer_init = explorer.replicate()
msm.translate(explorer_init, [-1.0, 0.0, 0.0]*unit.nanometers)


# ## Random shifts in cartesian coordinates

# In[6]:


msm.view([explorer_init, explorer])


# In[7]:


explorer.move.cartesian_shifts.show_parameters()


# In[8]:


explorer.move.cartesian_shifts.run()


# In[9]:


explorer.move.cartesian_shifts.atoms_moved


# In[10]:


msm.info(explorer, target='atom', selection=explorer.move.cartesian_shifts.atoms_moved)


# In[11]:


explorer.move.cartesian_shifts.shifts_moved


# In[12]:


msm.view([explorer_init, explorer])


# In[13]:


explorer.move.cartesian_shifts.set_parameters(selection='atom_name=="C"', mode_atoms='all',
                                              mode_steps='rmsd', step_size=0.3*unit.nanometers)


# In[14]:


explorer.move.cartesian_shifts.show_parameters()


# In[15]:


explorer.set_coordinates(initial_positions)


# In[16]:


explorer.move.cartesian_shifts()


# In[17]:


explorer.move.cartesian_shifts.atoms_moved


# In[18]:


explorer.move.cartesian_shifts.shifts_moved


# In[19]:


msm.view([explorer_init, explorer])


# ## Random shifts of dihedral angles

# In[20]:


explorer.set_coordinates(initial_positions)


# In[21]:


msm.view([explorer_init, explorer])


# In[22]:


msm.get_dihedral_angles(explorer, dihedral_angle='phi-psi')


# In[23]:


explorer.move.dihedral_shifts.show_parameters()


# In[24]:


explorer.move.dihedral_shifts()


# In[25]:


explorer.move.dihedral_shifts.quartets_moved


# In[26]:


explorer.move.dihedral_shifts.shifts_moved


# In[27]:


msm.get_dihedral_angles(explorer, dihedral_angle='phi-psi')


# In[28]:


msm.view([explorer_init, explorer])


# In[ ]:





# In[ ]:




