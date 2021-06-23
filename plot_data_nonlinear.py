import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import ZCBF_module
import two_dof_module
import yaml
import sys
import my_plot_module as mpm
import nonlinear_transformation_module as ntm

# Choose data to plot
#data_file =  ['sim_num0_exp0.txt','sim_num2_exp2.txt',  'sim_num1_exp1.txt']
#yaml_file = ['ZCBF_control_exp0.yaml','ZCBF_control_exp2.yaml', 'ZCBF_control_exp1.yaml']

data_file =  ['sim_num6_exp5_nominal.txt', 'sim_num5_exp5_nonlinear.txt']
yaml_file = ['ZCBF_control_exp5_nominal.yaml', 'ZCBF_control_exp5_nonlinear.yaml']

# Define custom colors for plots
c_black = (0,0,0)
c_blue = (0, .447, .741)
c_sky_blue = (0.3010, 0.745, 0.933)
c_orange = (0.85, 0.325, 0.098)
c_yellow = (0.929, 0.694, 0.125)
c_green = (0.466, 0.674, 0.188)
c_grey = (0.25, 0.25, 0.25)
c_purple = (0.494, 0.184, 0.556)
c_light_blue = (0.3010, 0.7450, 0.9330)

plot_color = [c_orange, c_orange, c_green]
contour_color = [c_green, c_sky_blue]
line_width = 2.5
contour_width = 4.0
plot_counter = -1

#fig, axes = plt.subplots(3,2)
fontsize = 24
for ii in range(len(data_file)):

    # Load parameter and simulation files
    with open(yaml_file[ii]) as file:
        param_data = yaml.full_load(file)
    sim_data = np.loadtxt(data_file[ii])


    # Collect data
    t = sim_data[:,0]
    y = sim_data[:,1:5]
    u = sim_data[:,5:7]

    # Define system constraints
    q_max = param_data['q_max'] 
    q_min = param_data['q_min']
    v_max = param_data['v_max']
    v_min = param_data['v_min']
    u_max = param_data['u_max']
    u_min = param_data['u_min']
    T_sim = param_data['T_sim']

    # Define 2DOF model parameters
    m_list = param_data['m_list']# link masses (kg)
    l_list = param_data['l_list']# link lengths(m)        
    f_list = param_data['f_list']# damping friction term
    I_list = [(1.0/12.0)*m_list[ii]*l_list[ii]**2 for ii in range(len(m_list)) ]  # moments of inertia
    grav = param_data['grav'] # gravity constant

    # Consctruct 2DOF system dynamics class
    two_dof_sys = two_dof_module.TwoDOFSystem(m_list,l_list,f_list,I_list,grav,param_data['gains'],param_data['ref'] ,param_data['nom_control'],kc=param_data['kc'], kg=param_data['kg'], km2 = param_data['km2'], c3 = param_data['c3'], c1 = param_data['c1'], c5= param_data['c5']) #

    # Define nonlinear transform for nonlinear constraints
    try:
        ctype0 = param_data['nonlinear_params']['constraint0']['type']
        ctype1 = param_data['nonlinear_params']['constraint1']['type']
        ctype0_params = param_data['nonlinear_params']['constraint0']['params']
        ctype1_params = param_data['nonlinear_params']['constraint1']['params']
        n_transform = ntm.NonlinearTransformation(two_dof_sys.n_joints, [ctype0, ctype1], [ctype0_params, ctype1_params])
        print('Nonlinear tranformation added')
    except:
        n_transform = ntm.NonlinearTransformation(two_dof_sys.n_joints)
        print('No nonlinear transformation')

    plot_class = mpm.PlotClass()
    contour_range_0 = np.linspace(-3.0, 3.0, num=100)
    contour_range_1 = np.linspace(-3.0, 3.0, num=100)

    # Plot position plot
    plot_counter += 1
    plt.figure(plot_counter)
    for jj in range(n_transform.n):
        plot_class.add_contour(plot_counter, n_transform.c[jj].eval, contour_range_0, contour_range_1, levels = [q_min[jj]], color=[contour_color[jj]], linestyle='dashed', linewidth=[contour_width])
        plot_class.add_contour(plot_counter, n_transform.c[jj].eval, contour_range_0, contour_range_1, levels = [q_max[jj]], color=[contour_color[jj]], linewidth=[contour_width])

    plot_class.plot_2D_phase(plot_counter, y[:,:2])
    #plt.plot(y[:,0], y,c=plot_color[ii],linewidth=line_width,label='q0')
    #plt.plot([0,T_sim], [q_max[0], q_max[0]], 'k--')
    #plt.plot([0,T_sim], [q_min[0], q_min[0]], 'k--')
    plt.tick_params(labelsize=fontsize)
    plt.xlim([0.0, 1.7])
    plt.ylim([1.3, 2.7])
    # Remove right and top box lines
    ax = plt.gca()
    right_side = ax.spines["right"]
    top_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side.set_visible(False)
    plt.xlabel(r'$\displaystyle \tilde{q}_0$', fontsize=fontsize)
    plt.ylabel(r'$\displaystyle \tilde{q}_1$', fontsize=fontsize)
    plt.tight_layout()

    # Plot v0
    plot_counter += 1
    plt.figure(plot_counter)
    plt.plot(t,y[:,2],c=plot_color[ii],linewidth=line_width,label='v0')
    plt.plot([0,T_sim], [v_max[0], v_max[0]], 'k--')
    plt.plot([0,T_sim], [v_min[0], v_min[0]], 'k--')
    plt.tick_params(labelsize=fontsize)
    plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
    plt.ylabel(r'$\displaystyle v_0$', fontsize=fontsize)
    #plt.legend()
    plt.tight_layout()
    # Remove right and top box lines
    ax = plt.gca()
    right_side = ax.spines["right"]
    top_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side.set_visible(False)

    # Plot v1
    plot_counter += 1
    plt.figure(plot_counter)
    plt.plot(t,y[:,3],c=plot_color[ii],linewidth=line_width,label='v1')
    plt.plot([0,T_sim], [v_max[1], v_max[1]], 'k--')
    plt.plot([0,T_sim], [v_min[1], v_min[1]], 'k--')
    plt.tick_params(labelsize=fontsize)
    plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
    plt.ylabel(r'$\displaystyle v_1$', fontsize=fontsize)
    plt.tight_layout()
    # Remove right and top box lines
    ax = plt.gca()
    right_side = ax.spines["right"]
    top_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side.set_visible(False)
    #plt.legend()

    # Plot u0
    plot_counter += 1
    plt.figure(plot_counter)
    plt.plot(t,u[:,0],c=plot_color[ii],linewidth=line_width,label='u0')
    plt.plot([0,T_sim], [u_max[0], u_max[0]], 'k--')
    plt.plot([0,T_sim], [u_min[0], u_min[0]], 'k--')
    plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
    plt.tick_params(labelsize=fontsize)
    plt.tight_layout()
    plt.ylabel(r'$\displaystyle u_0$', fontsize=fontsize)
    plt.tight_layout()
    # Remove right and top box lines
    ax = plt.gca()
    right_side = ax.spines["right"]
    top_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side.set_visible(False)
    #plt.legend()

    # Plot u1
    plot_counter += 1
    plt.figure(plot_counter)
    plt.plot(t,u[:,1],c=plot_color[ii],linewidth=line_width,label='u1')
    plt.plot([0,T_sim], [u_max[1], u_max[1]], 'k--')
    plt.plot([0,T_sim], [u_min[1], u_min[1]], 'k--')
    plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
    plt.ylabel(r'$\displaystyle u_1$', fontsize=fontsize)
    plt.tick_params(labelsize=fontsize)
    plt.tight_layout()
    # Remove right and top box lines
    ax = plt.gca()
    right_side = ax.spines["right"]
    top_side = ax.spines["top"]
    right_side.set_visible(False)
    top_side.set_visible(False)
    
plt.show()


