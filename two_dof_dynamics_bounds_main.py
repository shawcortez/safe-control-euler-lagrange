import numpy as np
import matplotlib.pyplot as plt
import ZCBF_module
import two_dof_module
import yaml
import sys

#========================================
# Define simlation scenario

sim_file = 'ZCBF_control_'+ sys.argv[1] + '.yaml' # 'ZCBF_cont_control_expi.yaml'
print('Evaluating bounds for experiment: '+ sim_file)

#========================================

# Load parameter file
with open(sim_file) as file:
	data = yaml.full_load(file)


# Define system constraints
q_max = data['q_max'] 
q_min = data['q_min']
v_max = data['v_max']
v_min = data['v_min']
u_max = data['u_max']
u_min = data['u_min']


# Define 2DOF model parameters
m_list = data['m_list']# link masses (kg)
l_list = data['l_list']# link lengths(m)		
f_list = data['f_list']# damping friction term
I_list = [(1.0/12.0)*m_list[ii]*l_list[ii]**2 for ii in range(len(m_list)) ]  # moments of inertia
grav= data['grav'] # gravity constant

# Consctruct 2DOF system dynamics class
two_dof_sys = two_dof_module.TwoDOFSystem(m_list,l_list,f_list,I_list,grav,data['gains'],data['ref'] ,data['nom_control']) #

# Compute dynamic bounds
q_min_delta = q_min - data['delta']*np.ones(2)
q_max_delta = q_max + data['delta']*np.ones(2)
two_dof_sys.u_max = u_max
two_dof_sys.compute_system_bounds(q_min_delta.tolist(), q_max_delta.tolist() , v_min[:], v_max[:],dq_= data['dq_bounds'], dv_ = data['dv_bounds'])