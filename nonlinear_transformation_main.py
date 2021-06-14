import numpy as np
import matplotlib.pyplot as plt
import ZCBF_module
import two_dof_module
from scipy.integrate import odeint
import nonlinear_transformation_module as ntm
import my_plot_module as mpm

# Define system constraints
q_max = [0.95, 0.5]
q_min = [0.0, -0.5] 
v_max = [1.5, 1.5]
v_min = [-1.5, -1.5]
u_max = [10.0, 10.0]
u_min = [-10.0, -10.0]

# Define 2DOF model parameters
m_list = [1.0, 1.0]     # link masses (kg)
l_list = [1.0, 1.0]     # link lengths(m)       
f_list = [0.5, 0.5] #damping friction term
I_list = [(1/12)*m_list[0]*l_list[0]**2, (1/12)*m_list[1]*l_list[1]**2 ] # moments of inertia
grav= 0.0;  # gravity constant


# Define ZCBF parameters
gamma = 1
nu = .1
delta = 0.2
eta_bar = 0.01

# Consctruct 2DOF system dynamics class
gains = { 'Kp': 1.0, 'Kv': 0.1}
two_dof_sys = two_dof_module.TwoDOFSystem(m_list,l_list,f_list,I_list,grav,gains,None,None)

# Construct nonlinear tranformation
params_ellipsoid = {}
params_ellipsoid['a'] = 1.0
params_ellipsoid['P']= [1.0, 0.0 , 0.0, -1.0] # nxn matrix flattened into list
params_ellipsoid['qr']= [-1.0, 0.0]  # center of ellipsoid
params_planar = {}
params_planar['a'] = 0.0
params_planar['P']= [0.0, 0.0 , 0.0, 0.0] # nxn matrix flattened into list
params_planar['qr']= [0.1, 1.0]  # center of ellipsoid
n_transform = ntm.NonlinearTransformation(two_dof_sys, ['ellipsoid', 'planar'], [params_ellipsoid, params_planar])

# Compute Gradient
print('Gradient: ')
print(n_transform.eval_gradient(np.array([0.0, 0.0])))
print('end')

print('Invert gradient: ')
print(np.linalg.inv(n_transform.eval_gradient(np.array([0.0, 0.0]))))

#Compute Hessian
#print(n_transform.eval_hessian(np.array([0.1, 0.1])))

#Multiply Hessian
#v_test = np.array([0.1, 0.1])
#print(v_test.dot( n_transform.eval_hessian(np.array([0.1, 0.1])).dot(v_test)   ))


plot_class = mpm.PlotClass()
contour_range = np.linspace(-3.0, 3.0, num=100)



fig_ct = 1

for ii in range(n_transform.n):
    plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_min[ii]], color='red')
    plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_max[ii]], color='orange')
    #plot_class.add_fill(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, level = [0.01, 0.05])



plt.show()