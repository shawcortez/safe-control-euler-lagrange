import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import ZCBF_module
import two_dof_module
import yaml
import sys

# Choose data to plot
#data_file =  ['sim_num0_exp0.txt','sim_num2_exp2.txt',  'sim_num1_exp1.txt']
#yaml_file = ['ZCBF_control_exp0.yaml','ZCBF_control_exp2.yaml', 'ZCBF_control_exp1.yaml']

data_file =  ['sim_num0_exp0.txt', 'sim_num3_exp2_discrete.txt']
yaml_file = ['ZCBF_control_exp0.yaml', 'ZCBF_control_exp2_discrete.yaml']

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

plot_color = [c_orange, c_sky_blue, c_green]
line_width = 2.5

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

	plot_counter = 0

	# Plot q0
	plt.figure(plot_counter)
	plt.plot(t,y[:,0],c=plot_color[ii],linewidth=line_width,label='q0')
	plt.plot([0,T_sim], [q_max[0], q_max[0]], 'k--')
	plt.plot([0,T_sim], [q_min[0], q_min[0]], 'k--')
	plt.tick_params(labelsize=fontsize)
	#plt.title(r'$\displaystyle q_0$', fontsize=fontsize)
	#plt.xlabel('time')
	#plt.ylabel('q0(t)')
	#plt.legend()


	# Plot q1
	plot_counter += 1
	plt.figure(plot_counter)
	plt.plot(t,y[:,1],c=plot_color[ii],linewidth=line_width,label='q1')
	plt.plot([0,T_sim], [q_max[1], q_max[1]], 'k--')
	plt.plot([0,T_sim], [q_min[1], q_min[1]], 'k--')
	plt.tick_params(labelsize=fontsize)
	#plt.xlabel('time')
	#plt.ylabel('q1(t)')
	#plt.legend()

	# Plot v0
	plot_counter += 1
	plt.figure(plot_counter)
	plt.plot(t,y[:,2],c=plot_color[ii],linewidth=line_width,label='v0')
	plt.plot([0,T_sim], [v_max[0], v_max[0]], 'k--')
	plt.plot([0,T_sim], [v_min[0], v_min[0]], 'k--')
	plt.tick_params(labelsize=fontsize)
	#plt.xlabel('time')
	#plt.ylabel('v0(t)')
	#plt.legend()

	# Plot v1
	plot_counter += 1
	plt.figure(plot_counter)
	plt.plot(t,y[:,3],c=plot_color[ii],linewidth=line_width,label='v1')
	plt.plot([0,T_sim], [v_max[1], v_max[1]], 'k--')
	plt.plot([0,T_sim], [v_min[1], v_min[1]], 'k--')
	plt.tick_params(labelsize=fontsize)
	#plt.xlabel('time')
	#plt.ylabel('v1(t)')
	#plt.legend()

	# Plot u0
	plot_counter += 1
	plt.figure(plot_counter)
	plt.plot(t,u[:,0],c=plot_color[ii],linewidth=line_width,label='u0')
	plt.plot([0,T_sim], [u_max[0], u_max[0]], 'k--')
	plt.plot([0,T_sim], [u_min[0], u_min[0]], 'k--')
	#plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
	plt.tick_params(labelsize=fontsize)
	#plt.ylabel('u(t)')
	#plt.legend()

	# Plot u1
	plot_counter += 1
	plt.figure(plot_counter)
	plt.plot(t,u[:,1],c=plot_color[ii],linewidth=line_width,label='u1')
	plt.plot([0,T_sim], [u_max[1], u_max[1]], 'k--')
	plt.plot([0,T_sim], [u_min[1], u_min[1]], 'k--')
	#plt.xlabel(r'$\displaystyle t$', fontsize=fontsize)
	plt.tick_params(labelsize=fontsize)
	#plt.ylabel('u(t)')
	#plt.legend()

	# axes[0,0].plot(t,y[:,0],c=plot_color[ii],linewidth=2,label='q0')
	# axes[0,1].plot(t,y[:,1],c=plot_color[ii],linewidth=2,label='q1')
	# axes[1,0].plot(t,y[:,2],c=plot_color[ii],linewidth=2,label='v0')
	# axes[1,1].plot(t,y[:,3],c=plot_color[ii],linewidth=2,label='v1')
	# axes[2,0].plot(t,u[:,0],c=plot_color[ii],linewidth=2,label='u0')
	# axes[2,1].plot(t,u[:,1],c=plot_color[ii],linewidth=2,label='u1')
	

	# axes[2,0].set_xlabel(r't', fontsize=fontsize)
	# axes[2,1].set_xlabel(r't', fontsize=fontsize)
	

	#axes[2,1].legend()

# axes[0,0].plot([0,T_sim], [q_max[0], q_max[0]], 'k--')
# axes[0,0].plot([0,T_sim], [q_min[0], q_min[0]], 'k--')
# axes[0,0].set_ylabel(r'$\displaystyle q_0$', fontsize=fontsize)
# axes[0,1].plot([0,T_sim], [q_max[1], q_max[1]], 'k--')
# axes[0,1].plot([0,T_sim], [q_min[1], q_min[1]], 'k--')
# axes[0,1].set_ylabel(r'$\displaystyle q_1$', fontsize=fontsize)
# axes[1,0].plot([0,T_sim], [v_max[0], v_max[0]], 'k--')
# axes[1,0].plot([0,T_sim], [v_min[0], v_min[0]], 'k--')
# axes[1,0].set_ylabel(r'$\displaystyle v_0$', fontsize=fontsize)
# axes[1,1].plot([0,T_sim], [v_max[1], v_max[1]], 'k--')
# axes[1,1].plot([0,T_sim], [v_min[1], v_min[1]], 'k--')
# axes[1,1].set_ylabel(r'$\displaystyle v_1$', fontsize=fontsize)
# axes[2,0].plot([0,T_sim], [u_max[0], u_max[0]], 'k--')
# axes[2,0].plot([0,T_sim], [u_min[0], u_min[0]], 'k--')
# axes[2,0].set_ylabel(r'$\displaystyle u_0$', fontsize=fontsize)
# axes[2,1].plot([0,T_sim], [u_max[1], u_max[1]], 'k--')
# axes[2,1].plot([0,T_sim], [u_min[1], u_min[1]], 'k--')
# axes[2,1].set_ylabel(r'$\displaystyle u_1$', fontsize=fontsize)

plt.show()


