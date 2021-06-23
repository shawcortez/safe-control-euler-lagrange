import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import ZCBF_module
import two_dof_module
import yaml
import sys
import math

class PlotClass():
	def __init__(self, fontsize = 20, line_width = 1.5):

		# Define custom colors for plots
		self.c_black = (0,0,0)
		self.c_blue = (0, .447, .741)
		self.c_sky_blue = (0.3010, 0.745, 0.933)
		self.c_orange = (0.85, 0.325, 0.098)
		self.c_yellow = (0.929, 0.694, 0.125)
		self.c_green = (0.466, 0.674, 0.188)
		self.c_grey = (0.25, 0.25, 0.25)
		self.c_purple = (0.494, 0.184, 0.556)
		self.c_light_blue = (0.3010, 0.7450, 0.9330)

		self.c_string = [self.c_black, self.c_blue, self.c_sky_blue, self.c_orange, self.c_yellow, self.c_green, self.c_grey, self.c_purple, self.c_light_blue]


		self.c_neon_green = (0.0, 1.0, 0.0)
		self.plot_counter = 0

		self.line_width = line_width
		self.fontsize = fontsize


	def plot_states_control(self,t,q,v,u):
		'''Plot 3x2 subplots of q, v, u over time'''
		self.plot_counter += 1
		fig, axes = plt.subplots(3,2)

		axes[0,0].plot(t,q[:,0],'b',linewidth=self.line_width,label='q0')
		axes[0,0].set_ylabel(r'$\displaystyle q_0$', fontsize=self.fontsize)
		#axes[0,0].legend()
		axes[0,1].plot(t,q[:,1],'b',linewidth=self.line_width,label='q1')
		axes[0,1].set_ylabel(r'$\displaystyle q_1$', fontsize=self.fontsize)
		#axes[0,1].legend()
		axes[1,0].plot(t,v[:,0],'b',linewidth=self.line_width,label='v0')
		axes[1,0].set_ylabel(r'$\displaystyle v_0$', fontsize=self.fontsize)
		#axes[1,0].legend()
		axes[1,1].plot(t,v[:,1],'b',linewidth=self.line_width,label='v1')
		axes[1,1].set_ylabel(r'$\displaystyle v_1$', fontsize=self.fontsize)
		#axes[1,1].legend()
		axes[2,0].plot(t,u[:,0],'b',linewidth=self.line_width,label='u0')
		axes[2,0].set_ylabel(r'$\displaystyle u_0$', fontsize=self.fontsize)
		#axes[2,0].legend()
		axes[2,1].plot(t,u[:,1],'b',linewidth=self.line_width,label='v1')
		axes[2,1].set_ylabel(r'$\displaystyle u_1$', fontsize=self.fontsize)

		axes[2,0].set_xlabel(r't', fontsize=self.fontsize)
		axes[2,1].set_xlabel(r't', fontsize=self.fontsize)

		return self.plot_counter

	def plot_2D_phase(self,plot_count, q, label_str= ['q_0', 'q_1'],colors='black'):
		'''Plot q1 vs q2 trajectory'''
		
		plt.figure(plot_count)
		plt.plot(q[:,0],q[:,1],c= colors, linewidth= self.line_width, label='traj')
		#plt.plot(q[0,0],q[0,1],c= self.c_green, marker='*',linewidth= 10*self.line_width, label='start')
		#plt.plot(q[-1,0], q[-1,1], c= self.c_orange, marker='*', linewidth= 10*self.line_width, label= 'end')
		# plt.tick_params(labelsize=self.fontsize)
		# plt.legend()

		# plt.xlabel(r'$\displaystyle '+label_str[0]+'$',fontsize=self.fontsize)
		# plt.ylabel(r'$\displaystyle '+label_str[1]+'$',fontsize=self.fontsize)
		# plt.axis('equal')

		return

	def add_2D_phase(self, plt_count, q, color='black'):
		'''Plot q1 vs q2 in existing phase plot'''

		plt.figure(plt_count)

		plt.plot(q[:,0],q[:,1],c= color, linewidth= self.line_width, linestyle='dashdot')
		#plt.plot(q[0,0],q[0,1],c= self.c_green, marker='*',linewidth= 10*self.line_width)
		#plt.plot(q[-1,0], q[-1,1], c= self.c_orange, marker='*', linewidth= 10*self.line_width, label= 'end')
		#plt.tick_params(labelsize=self.fontsize)
		#plt.legend()

		#plt.xlabel(r'$\displaystyle '+label_str[0]+'$',fontsize=self.fontsize)
		#plt.ylabel(r'$\displaystyle '+label_str[1]+'$',fontsize=self.fontsize)
		

	def add_destab_points(self, plt_count, t, q, v, func):
		for ii in range(len(t)):
			qi = q[ii,:]
			vi = v[ii,:]
			ti = t[ii]

			# Determine region
			reg = func(qi,vi)
			if reg == 'region_1':
				plt.plot(qi[0],qi[1],linewidth=self.line_width,color='red',marker='o',markersize=2)
			


	def add_contour(self,plt_count, func,x_plot, y_plot, levels=[0.0], color='green', linestyle='solid', linewidth=[1.5]):
		'''Add contour of give func function to an existing plot for specified levels'''
		n = len(x_plot)
		m = len(y_plot)
		
		for kk in range(len(levels)):
			z_plot = np.zeros((n,m))
			for ii in range(n):
				for jj in range(m):
					var = np.array([x_plot[ii], y_plot[jj]])
					z_plot[ii,jj] = func(var)

			plt.figure(plt_count)
			plt.contour(x_plot, y_plot, z_plot.T, levels[kk], colors=color, linestyles=linestyle, linewidths= linewidth)

	def add_2Dgrid(self,plt_count, q_min, q_max, dq, color='black'):
		'''Add grid points defined by q_min, q_max with interval dq'''

		q0_vec = np.linspace(q_min[0],q_max[0], 1/dq[0])
		q1_vec = np.linspace(q_min[1],q_max[1], 1/dq[1])

		plt.figure(plt_count)

		for qi in q0_vec:
			for qj in q1_vec:
				plt.plot(qi,qj,color=color,marker='o',markersize=0.5)


	def add_fill(self, plt_count, func, x_plot, y_plot,level):
		'''Add fill depending on function sign'''
		n = len(x_plot)
		m = len(y_plot)
		z_plot = np.zeros((n,m))
		for ii in range(n):
			for jj in range(m):
				var = np.array([x_plot[ii], y_plot[jj]])
				z_plot[ii,jj] = func(var)

		plt.figure(plt_count)
		ax = plt.gca()
		#cs = plt.contourf(x_plot, y_plot, z_plot.T)
		#plt.colorbar(cs)
		cs2 = plt.contourf(x_plot, y_plot, z_plot.T,levels=level, extend='both')
		plt.colorbar(cs2)

	def add_2d_to_scalar_trajectory(self,plt_count, t,q,v,func,label='func(q(t),v(t))',color='black'):

		# Evaluate function over state trajectory
		h = []
		for ii in range(len(t)):
			qi = q[ii,:]
			vi = v[ii,:]
			h.append(func(qi,vi) )

		plt.figure(plt_count)
		plt.plot(t,h,linewidth=self.line_width,color=color)

		return
        
	def plot_scalar_trajectory(self,t,x,func,label='func(t)'):
		'''Plot a function of a state trajectory x over t'''

		self.plot_counter += 1

		# Evaluate function over state trajectory
		h = np.array([ func(xi) for xi in x ])

		plt.figure(self.plot_counter)
		plt.plot(t,h,linewidth=self.line_width)
		plt.tick_params(labelsize=self.fontsize)

		plt.xlabel(r'$\displaystyle t(0)$',fontsize=self.fontsize)
		plt.ylabel(r'$\displaystyle '+label+'$',fontsize=self.fontsize)

		return self.plot_counter

	def plot_2d_to_scalar_trajectory(self,t,q,v,func,label='func(q(t),v(t))',color='black'):

		self.plot_counter += 1

		# Evaluate function over state trajectory
		h = []
		for ii in range(len(t)):
			qi = q[ii,:]
			vi = v[ii,:]
			h.append(func(qi,vi) )

		plt.figure(self.plot_counter)
		plt.plot(t,h,linewidth=self.line_width,color=color)
		plt.tick_params(labelsize=self.fontsize)

		plt.xlabel(r'$\displaystyle t(0)$',fontsize=self.fontsize)
		plt.ylabel(r'$\displaystyle '+label+'$',fontsize=self.fontsize)

		return

	def plot_Lyapunov_with_regions(self,t,q,v,Lyapunov_func, stab_region_func,label):
		self.plot_counter += 1

		# Evaluate function over state trajectory
		
		plt.figure(self.plot_counter)
		for ii in range(len(t)):
			qi = q[ii,:]
			vi = v[ii,:]
			ti = t[ii]

			# Compute V
			Vi = Lyapunov_func(qi,vi) 

			# Determine region
			reg = stab_region_func(qi,vi)
			if reg == 'region_1':
				plt.plot(ti,Vi,linewidth=self.line_width,color='red',marker='o',markersize=2)
			if reg == 'region_2':
				plt.plot(ti,Vi,linewidth=self.line_width,color='green',marker='o',markersize=2)
			if reg == 'region_3':
				plt.plot(ti,Vi,linewidth=self.line_width,color='black',marker='o',markersize=2)
			if reg == 'region_4':
				plt.plot(ti,Vi,linewidth=self.line_width,color='orange',marker='o',markersize=2)

		
		
		plt.tick_params(labelsize=self.fontsize)

		plt.xlabel(r'$\displaystyle t(s)$',fontsize=self.fontsize)
		plt.ylabel(r'$\displaystyle '+label+'$',fontsize=self.fontsize)

		return