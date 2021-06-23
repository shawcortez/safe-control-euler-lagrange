import numpy as np
import matplotlib.pyplot as plt
import ZCBF_module
import two_dof_module
from scipy.integrate import odeint
import nonlinear_transformation_module as ntm
import my_plot_module as mpm
import yaml
import sys

#========================================
# Define simlation scenario

sim_file = 'ZCBF_control_'+ sys.argv[1] + '.yaml' # 'ZCBF_cont_control_expi.yaml'
print('Running experiment: '+ sim_file)

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
two_dof_sys = two_dof_module.TwoDOFSystem(m_list,l_list,f_list,I_list,grav,data['gains'],data['ref'] ,data['nom_control'],kc=data['kc'], kg=data['kg'], km2 = data['km2'], c3 = data['c3'], c1 = data['c1'], c5= data['c5']) #

# Define nonlinear transform for nonlinear constraints
try:
    ctype0 = data['nonlinear_params']['constraint0']['type']
    ctype1 = data['nonlinear_params']['constraint1']['type']
    ctype0_params = data['nonlinear_params']['constraint0']['params']
    ctype1_params = data['nonlinear_params']['constraint1']['params']
    n_transform = ntm.NonlinearTransformation(two_dof_sys.n_joints, [ctype0, ctype1], [ctype0_params, ctype1_params])
    print('Nonlinear tranformation added')
except:
    n_transform = ntm.NonlinearTransformation(two_dof_sys.n_joints)
    print('No nonlinear transformation')

# Construct ZCBF class
# Define ZCBF (initial) parameters
gamma = data['gamma'] 
nu = data['nu']
delta = data['delta'] 
eta_bar = data['eta_bar'] 
alpha = data['alpha']
beta = data['beta']
zcbf_el1 = ZCBF_module.ZCBF(alpha,beta,q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,two_dof_sys,n_transform)




# Compute hessian bound:
zcbf_el1.compute_nonlinear_fx_bound( data['nonlinear_params']['qt_min'], data['nonlinear_params']['qt_max'], dq=data['dq'])

# Compute dynamic bounds
#zcbf_el1.compute_system_bounds(data['nonlinear_params']['qt_min'], data['nonlinear_params']['qt_max'] , v_min[:], v_max[:],dq_= data['dq_bounds'], dv_ = data['dv_bounds'])

# Test hessian computations
# x = np.array([0.0, 2.0])
# print('x: '+ str(x))
# grad = n_transform.eval_gradient(x)
# print('Gradient: '+str(grad))
# H = n_transform.eval_hessian(x)
# print('Hessian: '+ str(H))
# v = np.array([1.0, 1.0])
# print('v: '+ str(v))
# print('Hv = ' + str( H.dot(v) ))
# print('vTHv = '+ str( v.dot(H.dot(v) )))

plot_class = mpm.PlotClass()
contour_range = np.linspace(-3.0, 3.0, num=100)



fig_ct = 1

for ii in range(n_transform.n):
    plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_min[ii]], color='red')
    plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_max[ii]], color='orange')
    #plot_class.add_fill(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, level = [0.01, 0.05])



plt.show()