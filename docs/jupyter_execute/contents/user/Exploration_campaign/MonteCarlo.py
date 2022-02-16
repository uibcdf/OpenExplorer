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
import matplotlib.pyplot as plt


# # MonteCarlo
# 
# Metropolis, N., Rosenbluth, A. W., Rosenbluth, M. N., Teller, A. H., & Teller, E. (1953). Equation of state calculations by fast computing machines. The journal of chemical physics, 21(6), 1087-1092.
# 
# Metropolis, N., & Ulam, S. (1949). The monte carlo method. Journal of the American statistical association, 44(247), 335-341.
# 
# ## Goal
# 
# - Global optimization
# - Free energy exploration
# - Equilibrium magnitudes and distribution
# - Thermodinamic magnitudes

# ## Exploration Strategy

# In[3]:


#modeller = msm.convert('alanine_dipeptide.pdb', to_form='openmm.Modeller')
modeller = msm.convert('alanine_tetrapeptide.pdb', to_form='openmm.Modeller')

topology = modeller.topology
positions = modeller.positions

forcefield = app.ForceField('amber10.xml', 'amber10_obc.xml')

system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, #constraints=app.HBonds,
                                implicitSolvent=app.OBC2, soluteDielectric=1.0,
                                 solventDielectric=78.5)


# In[ ]:


explorer = oe.Explorer(topology, system, platform='CUDA')
explorer.set_coordinates(positions)


# In[ ]:


quartets = msm.covalent_dihedral_quartets(explorer, dihedral_angle='all')
print(msm.get_dihedral_angles(explorer, quartets=quartets))


# In[ ]:


print(explorer.get_potential_energy())


# - random move
# - accept

# In[ ]:


explorer.move.dihedral_shifts.set_parameters(dihedral_angle='all', mode_angles='random',
                                             n_random_angles=1, stepsize=1.5*unit.degrees,
                                             mode_steps='random')


# In[ ]:


explorer.acceptance.metropolis_hastings.set_parameters(temperature=500.0*unit.kelvin)


# In[ ]:


n_montecarlo_moves = 250000

traj = [] * unit.degrees
global_minimum_potential_energy = np.inf * unit.kilojoules_per_mole

coordinates = explorer.get_coordinates()
potential_energy = explorer.get_potential_energy()
global_minimum_potential_energy = potential_energy
global_minimum_coordinates = coordinates

angles = msm.get_dihedral_angles(explorer, quartets=quartets)
traj.append(angles)

previous_coordinates = coordinates
previous_potential_energy = potential_energy

for _ in tqdm(range(n_montecarlo_moves)):
    
    explorer.move.dihedral_shifts()
    coordinates = explorer.get_coordinates()
    potential_energy = explorer.get_potential_energy()
    explorer.acceptance.metropolis_hastings(previous_coordinates=previous_coordinates,
                                           previous_potential_energy=previous_potential_energy,
                                           coordinates=coordinates,
                                           potential_energy=potential_energy)
    
    if explorer.acceptance.metropolis_hastings.accepted:
        previous_coordinates = coordinates
        previous_potential_energy = potential_energy
        angles = msm.get_dihedral_angles(explorer, quartets=quartets)
        if potential_energy < global_minimum_potential_energy:
            global_minimum_potential_energy = potential_energy
            global_minimum_coordinates = coordinates

    traj.append(angles)


# In[ ]:


explorer.acceptance.metropolis_hastings.n_tries


# In[ ]:


explorer.acceptance.metropolis_hastings.n_accepted


# In[ ]:


global_minimum_potential_energy


# In[ ]:


traj._value = np.array(traj._value)


# In[ ]:


print(traj.shape)


# In[ ]:


plt.scatter(traj[:,0,0], traj[:,0,1])
plt.show()


# In[ ]:


import seaborn as sns

ax = sns.kdeplot(traj[:,0,0], traj[:,0,1], shade=True)
ax.set_xlim(-180.0,180.0)
ax.set_ylim(-180.0,180.0)
ax.set_xticks([-180.0, -90.0, 0.0, 90.0, 180.0])
ax.set_yticks([-180.0, -90.0, 0.0, 90.0, 180.0])
ax.set_xlabel('$\phi_1$')
ax.set_ylabel('$\psi_1$')
ax.set_aspect('equal')


# In[ ]:




