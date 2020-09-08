import numpy as np
import matplotlib.pyplot as plt
import ZCBF_module
import two_dof_module
from scipy.integrate import odeint

# Define system constraints
q_max = [1.0, 1.0]
q_min = [-1.0, -1.0] 
v_max = [1.5, 1.5]
v_min = [-1.5, -1.5]
u_max = [10.0, 10.0]
u_min = [-10.0, -10.0]

# Define 2DOF model parameters
m_list = [1.0, 1.0]  	# link masses (kg)
l_list = [1.0, 1.0]		# link lengths(m)		
f_list = [0.5, 0.5]	#damping friction term
I_list = [(1/12)*m_list[0]*l_list[0]**2, (1/12)*m_list[1]*l_list[1]**2 ] # moments of inertia
grav= 0.0;	# gravity constant


# Define ZCBF parameters
gamma = 1
nu = .1
delta = 0.2
eta_bar = 0.01

# Consctruct 2DOF system dynamics class
gains = { 'Kp': 1.0, 'Kv': 0.1}
two_dof_sys = two_dof_module.TwoDOFSystem(m_list,l_list,f_list,I_list,grav,gains)

# Construct ZCBF class with \alpha(h) = arctan(h), \beta(h) = h^3
zcbf_el1 = ZCBF_module.ZCBF('atan','cubic',q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,two_dof_sys)
#zcbf_el2 = ZCBF_module.ZCBF('linear','cubic',q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,two_dof_sys)

# Plot H sets
zcbf_el1.plot_H()
#zcbf_el2.plot_H()
