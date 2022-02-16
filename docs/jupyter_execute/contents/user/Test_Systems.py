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


# # Test systems

# ## Alanine dipeptide

# In[3]:


msm.build_peptide('aminoacids3:AceAlaNme', to_form='alanine_dipeptide.pdb', verbose=False)


# In[8]:


## Alanine tetrapeptide


# In[9]:


molecular_system = msm.build_peptide('aminoacids3:AceAlaAlaAlaNme', to_form='alanine_tetrapeptide.pdb', verbose=False)


# In[10]:


## Met-Enkephalin


# In[11]:


molecular_system = msm.build_peptide('aminoacids3:TyrGlyGlyPheMet', to_form='molsysmt.MolSys', verbose=False)


# In[12]:


molecular_system = msm.terminal_capping(molecular_system)


# In[13]:


msm.convert(molecular_system, to_form='metenkephalin.pdb')


# In[ ]:




