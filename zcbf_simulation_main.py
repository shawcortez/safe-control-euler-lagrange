import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import ZCBF_module
import two_dof_module
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
try:
    zcbf_el1 = ZCBF_module.ZCBF(alpha,beta,q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,two_dof_sys,n_transform, kh=data['nonlinear_params']['kh'])
except:
    zcbf_el1 = ZCBF_module.ZCBF(alpha,beta,q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,two_dof_sys,n_transform)


# Compute safe-by-design ZCBF parameters
T_sample = data['T_sample']
if data['compute_zcbf']:
    if n_transform.id_bool:
	   zcbf_el1.compute_zcbf_parameters(data['sim_control'], q_min, q_max, dq=data['dq'], T_sample=T_sample, eps_frac=data['eps_frac'], nu_frac= data['nu_frac'])
    else: 
        zcbf_el1.compute_zcbf_parameters(data['sim_control'], data['nonlinear_params']['qt_min'], data['nonlinear_params']['qt_max'], dq=data['dq'], T_sample=T_sample, eps_frac=data['eps_frac'], nu_frac= data['nu_frac'],nt_id_bool=False)
#raise(SystemExit)

# Setup ode solver
q0 = data['q0'] 
v0 = data['v0']
x0 = np.array(q0+v0)
T_sim = data['T_sim']
N_sim = data['N_sim']
t_sim = np.linspace(0.0,T_sim,N_sim)

#--------------------------------
# Solve odeint and run simulation

if data['sim_control'] == 'cont_time':
	# Continous control 
	t,y,u = zcbf_el1.run_simulation(t_sim,x0,data['cont_type'],rtol=1e-6,atol=1e-7,hmax=0.01,cont=True,T=None)

if data['sim_control'] == 'discrete_time':
	# Sampled-data control
	t,y,u = zcbf_el1.run_simulation(t_sim,x0,'ZCBF_control_sampled',rtol=1e-6,atol=1e-7,hmax=0.01,cont=False,T=T_sample)


 

fig, axes = plt.subplots(3,2)
fontsize = 15
axes[0,0].plot(t,y[:,0],'b',linewidth=2,label='q0')
axes[0,0].plot([0,T_sim], [q_max[0], q_max[0]], 'k--')
axes[0,0].plot([0,T_sim], [q_min[0], q_min[0]], 'k--')
axes[0,0].set_ylabel(r'$\displaystyle q_0$', fontsize=fontsize)
#axes[0,0].legend()
axes[0,1].plot(t,y[:,1],'b',linewidth=2,label='q1')
axes[0,1].plot([0,T_sim], [q_max[1], q_max[1]], 'k--')
axes[0,1].plot([0,T_sim], [q_min[1], q_min[1]], 'k--')
axes[0,1].set_ylabel(r'$\displaystyle q_1$', fontsize=fontsize)
#axes[0,1].legend()
axes[1,0].plot(t,y[:,2],'b',linewidth=2,label='v0')
axes[1,0].plot([0,T_sim], [v_max[0], v_max[0]], 'k--')
axes[1,0].plot([0,T_sim], [v_min[0], v_min[0]], 'k--')
axes[1,0].set_ylabel(r'$\displaystyle v_0$', fontsize=fontsize)
#axes[1,0].legend()
axes[1,1].plot(t,y[:,3],'b',linewidth=2,label='v1')
axes[1,1].plot([0,T_sim], [v_max[1], v_max[1]], 'k--')
axes[1,1].plot([0,T_sim], [v_min[1], v_min[1]], 'k--')
axes[1,1].set_ylabel(r'$\displaystyle v_1$', fontsize=fontsize)
#axes[1,1].legend()
axes[2,0].plot(t,u[:,0],'b',linewidth=2,label='u0')
axes[2,0].plot([0,T_sim], [u_max[0], u_max[0]], 'k--')
axes[2,0].plot([0,T_sim], [u_min[0], u_min[0]], 'k--')
axes[2,0].set_ylabel(r'$\displaystyle u_0$', fontsize=fontsize)
#axes[2,0].legend()
axes[2,1].plot(t,u[:,1],'b',linewidth=2,label='v1')
axes[2,1].plot([0,T_sim], [u_max[1], u_max[1]], 'k--')
axes[2,1].plot([0,T_sim], [u_min[1], u_min[1]], 'k--')
axes[2,1].set_ylabel(r'$\displaystyle u_1$', fontsize=fontsize)

axes[2,0].set_xlabel(r't', fontsize=fontsize)
axes[2,1].set_xlabel(r't', fontsize=fontsize)
#plt.tick_params(labelsize=15)

#axes[2,1].legend()

plot_class = mpm.PlotClass()
contour_range = np.linspace(-3.0, 3.0, num=100)



fig_ct = 2
if n_transform == None:
    plot_class.add_contour(fig_ct, ntm.eval_id, contour_range, contour_range, levels = [q_min[ii]], color='red')
    plot_class.add_contour(fig_ct, ntm.eval_id, contour_range, contour_range, levels = [q_max[ii]], color='orange')

else:
    for ii in range(zcbf_el1.n):
        plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_min[ii]], color='red')
        plot_class.add_contour(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, levels = [q_max[ii]], color='orange')
        #plot_class.add_fill(fig_ct, n_transform.c[ii].eval, contour_range, contour_range, level = [0.01, 0.05])
plot_class.add_2D_phase(fig_ct, y[:,0:2], color='black')


plt.show()

# Save data?
save_me = raw_input("Save data? (yes/no): ")
if save_me == 'yes':
	data_title = raw_input('Simulation number: ')
	sim_data = np.column_stack((t,y,u))
	np.savetxt('sim_num' + data_title + '_' + sys.argv[1] + '.txt', sim_data, header='Experiment: ' + str(sim_file) + ' \n  Time (1), State (' + str(2*zcbf_el1.n) + '), Input (' + str(zcbf_el1.n) +')')