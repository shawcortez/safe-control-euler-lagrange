import numpy as np
import math
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import two_dof_module
import cvxopt
from scipy.integrate import odeint


class ExtendedClassK():
    def __init__(self,type):
        
        # Define type of class-K
        if type in ["linear","cubic","atan"]:
            self.type = type
        else:
            raise TypeError('Extended Class K not recognized. (ExtendedClassK class)') 

            
    def eval(self,h):
        '''Evaluate the class-K for some input h'''
        if self.type == 'linear':
            return h
        if self.type == 'cubic':
            return h**3
        if self.type == 'atan':
            return np.arctan(h)
        
    def diff(self,h):
        if self.type == 'linear':
            return 1.0
        if self.type == 'cubic':
            return 3.0*h**2
        if self.type == 'atan':
            return 1.0/(h**2 + 1.0)

    def eval_L(self,hmax=None):
        if self.type == 'linear':
            return 1.0

        if self.type == 'cubic':
            return 3.0*hmax*2
        
        if self.type == 'atan':
            return 1.0

    def eval_zeta_delta(self,delta,gamma,q_max=None,q_min=None):
        if self.type == 'linear':
            return -gamma*delta

        if self.type == 'cubic':
            zeta = []
            for ii in range(len(q_max)):
                zeta.append( gamma*( self.eval( q_max[ii] - q_min[ii] + delta ) - self.eval(q_max[ii] - q_min[ii] + 2.0*delta)  ) )

            return min(zeta)

        if self.type == 'atan':
            return -gamma*2*self.eval(0.5*delta)

    def eval_rho_delta(self,delta,gamma,q_max,q_min,rho=None):
        if self.type == 'linear':
            return 0.5*gamma*min(np.array(q_max) - np.array(q_min))

        if self.type == 'atan':
            rho_qi = []

            # Min value of rho (for atan) occurs at boundary points where q_i = q_max_i + delta, q_i = q_min_i - delta
            for ii in range(len(q_max)):
                rho_upper_i = rho(ii,q_max[ii] + delta)
                rho_lower_i = rho(ii,q_min[ii] - delta)
                rho_qi.append(min(rho_upper_i, rho_lower_i))
            return min(rho_qi) 
            
class ZCBF():
    def __init__(self,alpha_type,beta_type,q_min,q_max,v_min,v_max,u_min,u_max,gamma,nu,delta,eta_bar,sys_dyn=None):
        '''Initialize ZCBF class with chosen \alpha, \beta functions, state/input constraints, and robustness margins'''
        
        # Store ZCBF parameters
        self.alpha = ExtendedClassK(alpha_type)
        self.beta = ExtendedClassK(beta_type)
        self.gamma = gamma
        self.nu = nu
        self.delta = delta
        self.eta = eta_bar
        
        # Store state/input constraint parameters
        self.q_min = q_min
        self.q_max = q_max
        self.v_min = v_min
        self.v_max = v_max
        self.u_min = u_min
        self.u_max = u_max
        self.n = len(q_min) # number of states in the system. Note q_min,q_max..., u_max are all of the same size.
        self.sys_dyn = sys_dyn

        #Check that q_min...u_max are appropriately defined
        if (self.n == len(q_max) and self.n == len(v_min) and self.n== len(v_max) and
                self.n == len(u_min) and self.n == len(u_max)):
            pass
        else:
            raise TypeError('q_min...,u_max vectors not of same length. (ZCBF class)') 

        # Define custom colors for plots
        self.c_black = (0,0,0)
        self.c_blue = (0, .447, .741)
        self.c_orange = (0.85, 0.325, 0.098)
        self.c_yellow = (0.929, 0.694, 0.125)
        self.c_green = (0.466, 0.674, 0.188)
        self.c_grey = (0.25, 0.25, 0.25)
        self.c_purple = (0.494, 0.184, 0.556)
        self.c_light_blue = (0.3010, 0.7450, 0.9330)
        
    def h_bar(self, ii, q_i):
        '''Define bar{h}_i (delta = 0)'''
        return self.q_max[ii] - q_i
        
    def h_ubar(self,ii,q_i):
        '''Define ubar{h}_i (delta = 0)'''
        return q_i - self.q_min[ii]
    
    def h_bar_delta(self,ii,q_i):
        '''Define bar{h}_i (delta \neq 0)'''
        return self.h_bar(ii,q_i) + self.delta
    
    def h_ubar_delta(self,ii,q_i):
        '''Define ubar{h}_i (delta \neq 0)'''
        return self.h_ubar(ii,q_i) + self.delta
    
    def b_bar(self,ii,q_i,v_i):
        '''Define bar{b}_i (delta = 0)'''
        return -v_i + self.gamma*self.alpha.eval(self.h_bar(ii,q_i))
    
    def b_ubar(self,ii,q_i,v_i):
        '''Define ubar{b}_i (delta = 0)'''
        return v_i + self.gamma*self.alpha.eval(self.h_ubar(ii,q_i))
    
    def b_bar_delta(self,ii,q_i,v_i):
        '''Define bar{b}_i (delta \neq 0)'''
        return -v_i + self.gamma*self.alpha.eval(self.h_bar_delta(ii,q_i))
    
    def b_ubar_delta(self,ii,q_i,v_i):
        '''Define ubar{b}_i (delta \neq 0)'''
        return v_i + self.gamma*self.alpha.eval(self.h_ubar_delta(ii,q_i))

    def v_upper(self,ii,q_i):
        '''Compute upper bound of v_i with delta =0'''
        return self.b_bar(ii,q_i,0)

    def v_upper_delta(self,ii,q_i):
        '''Compute upper bound of v_i with delta'''
        return self.b_bar_delta(ii,q_i,0)

    def v_lower(self,ii,q_i):
        '''Compute lower bound of v_i with delta =0'''
        return -self.b_ubar(ii,q_i,0)

    def v_lower_delta(self,ii,q_i):
        '''Compute lower bound of v_i with delta'''
        return -self.b_ubar_delta(ii,q_i,0)

    def rho(self,ii,q_i):
        '''Compute rho which satisfies b_bar + b_ubar = 2 rho in the safe set H'''
        return 0.5*self.gamma*( self.alpha.eval(self.h_bar(ii,q_i)) + self.alpha.eval(self.h_ubar(ii,q_i)) )

    def v_rho(self,ii,q_i):
        '''Compute the velocity when b_bar = b_ubar = rho '''
        return 0.5*self.gamma*( self.alpha.eval(self.h_bar(ii,q_i)) - self.alpha.eval(self.h_ubar(ii,q_i)) )

    def plot_H(self,index=0,n_length=100,delta_terms=True):
        '''Create H plot i.e the safe set defined by the barrier functions'''

        # Setup definitions for upper/lower bounds of H, H^delta, and rho manifold
        ax = plt.subplot(111)    
        H_upper_bound = []
        H_lower_bound = []
        H_upper_bound_delta = []
        H_lower_bound_delta = []
        rho_man = []
        q_i = np.linspace(self.q_min[index]-self.delta, self.q_max[index]+self.delta,n_length)
        for q in q_i:
            # Define upper and lower boundaries of H
            H_upper_bound.append(self.v_upper(index,q))
            H_lower_bound.append(self.v_lower(index,q))

            # Define manifold where b_bar = b_ubar = rho
            rho_man.append(self.v_rho(index,q))

            # Define upper and lower boundaries of H^delta
            H_upper_bound_delta.append(self.v_upper_delta(index,q))
            H_lower_bound_delta.append(self.v_lower_delta(index,q))

        # Define left and right boundaries of H
        q_left_bound = self.q_min[index]*np.ones(n_length)
        H_left_bound = np.linspace(min(self.v_min[index],min(H_lower_bound_delta)), max(self.v_max[index],max(H_upper_bound_delta)),n_length)
        q_right_bound = self.q_max[index]*np.ones(n_length)
        H_right_bound = H_left_bound

        # Define left and right boundaries of H^delta
        if delta_terms == True:
            q_left_bound_delta = (self.q_min[index]-self.delta)*np.ones(n_length)
            H_left_bound_delta = H_left_bound
            q_right_bound_delta = (self.q_max[index]+self.delta)*np.ones(n_length)
            H_right_bound_delta = H_right_bound

        # Plot boundaries
        plt.plot(q_i,H_upper_bound,'k--',linewidth=2.0)
        plt.plot(q_i,H_lower_bound,'k--',linewidth=2.0)
        plt.plot(q_left_bound,H_left_bound,'k--',linewidth=2.0)
        plt.plot(q_right_bound,H_right_bound,'k--',linewidth=2.0)
        if delta_terms == True:
            plt.plot(q_i,H_upper_bound_delta,'k',linewidth=2.0)
            plt.plot(q_i,H_lower_bound_delta,'k',linewidth=2.0)
            plt.plot(q_left_bound_delta,H_left_bound_delta,'k',linewidth=2.0)
            plt.plot(q_right_bound_delta,H_right_bound_delta,'k',linewidth=2.0)


        # Fill in area defining regions of H^delta
        v_axis = np.zeros(n_length)
        I1_lb = [max(y,0) for y in rho_man]
        I2_ub = [min(y,0) for y in H_upper_bound ]
        I3_lb = [max(y,0) for y in H_lower_bound]
        I4_ub = [min(y,0) for y in rho_man]

        # Region I
        ax.fill_between(q_i,I1_lb,H_upper_bound, where= H_upper_bound >= v_axis,facecolor = self.c_grey )
        # Region II
        ax.fill_between(q_i,rho_man,I2_ub, where= rho_man <= v_axis,facecolor = self.c_blue )
        # Region III
        ax.fill_between(q_i,I3_lb,rho_man, where= rho_man >= v_axis,facecolor = self.c_light_blue )
        # Region IV
        ax.fill_between(q_i,H_lower_bound,I4_ub, where= H_lower_bound <= v_axis,facecolor = self.c_orange )
        # Region V
        ax.fill_between(q_i,H_upper_bound_delta,H_upper_bound, facecolor = self.c_green )
        # Region VI
        ax.fill_between(q_i,H_lower_bound_delta,H_lower_bound, facecolor = self.c_yellow )
        # Region VII
        plt.plot(q_i[:n_length//2],rho_man[:n_length//2],'r-.',linewidth=2.0)
        # Region VIII
        plt.plot(q_i[n_length//2:],rho_man[n_length//2:],'w-.',linewidth=2.0)

        # Label plot components
        #plt.text(0, -1.25, r'$\displaystyle \mathcal{C}_i$', fontsize=16)
        

        # Define labels and show plot
        ax.set_xlabel(r'$\displaystyle q_i$', fontsize=24)
        ax.set_ylabel(r'$\displaystyle v_i$', fontsize=24)
        plt.tick_params(labelsize=15)
        plt.grid(True)
        plt.show()

    def determine_region_H(self,q,v):
        '''Determine what region of the safe set the current state is in'''
        for ii in range(self.n):
            qi = q[ii]
            vi = v[ii]
            b_bar = self.b_bar(ii,qi,vi)
            b_ubar = self.b_ubar(ii,qi,vi)
            rho = self.rho(ii,qi)
            if (b_bar >= 0.0) and (b_bar < rho) and (vi >= 0.0):
                return 1
            if (b_bar >= 0.0) and (b_bar < rho) and (vi <= 0.0):
                return 2
            if (b_ubar >= 0.0) and (b_ubar < rho) and (vi >= 0.0):
                return 3
            if (b_ubar >= 0.0) and (b_ubar < rho) and (vi <= 0.0):
                return 4
            if (b_bar < 0.0):
                return 5
            if (b_ubar < 0.0):
                return 6
            if (b_ubar == rho) and (vi >= 0.0):
                return 7
            if (b_ubar == rho) and (vi <= 0.0):
                return 8


    def compute_ubar_control(self,q,v):
        '''Compute available control law ubar that should always exist'''
        # Compute system dynamics
        M = self.sys_dyn.M_eval(q)
        C = self.sys_dyn.C_eval(q,v)
        g = self.sys_dyn.g_eval(q)

        mu = np.zeros(self.n)
        chi = np.zeros(self.n)
        psi = np.zeros(self.n)
        for ii in range(self.n):
            qi = q[ii]
            vi = v[ii]

            #determine region
            region = self.determine_region_H(q,v)

            # compute mu
            if (region == 1) or (region == 5) or (region == 7):
                mu[ii] = -self.gamma*self.alpha.diff(self.h_bar(ii,qi))*vi
            if (region == 2) or (region == 3):
                mu[ii] = 0.0
            if (region == 4) or (region == 6) or (region == 8):
                mu[ii] = -self.gamma*self.alpha.diff(self.h_ubar(ii,qi))*vi

            # compute chi
            if (region == 5):
                chi[ii] = self.nu*self.beta.eval(self.b_bar(ii,qi,vi))
            if (region == 6):
                chi[ii] = -self.nu*self.beta.eval(self.b_ubar(ii,qi,vi))
            # otherwise chi = 0.0

            # compute psi
            if (region == 1) or (region == 2) or (region == 5):
                psi[ii] = -self.eta
            if (region == 3) or (region == 4) or (region == 6):
                psi[ii] = self.eta
            #otherwise psi = 0.0

        return np.matmul(M, mu + chi + psi) + np.matmul(C,v)+ np.matmul(self.sys_dyn.F, v) - g





    def compute_zcbf_control(self,q,v,f_x,g_x,t):
        ''' Compute u*, the zcbf-based control law, given the 
        system dynamics f_x, g_x and nominal control law u_nom'''


        # Construct Lie derivatives of b w.r.t system dynamics
        S = np.concatenate((-np.eye(self.n),np.eye(self.n)))
        Lf_b1 = S.dot(f_x)
        Lg_b = S.dot(g_x)

        # Compute nominal control law
        u_nom = self.sys_dyn.nominal_control(q,v,t)

        # Construct zcbf based components of the QP
        zcbf_p = np.zeros(2*self.n)
        Lf_b2 = np.zeros(2*self.n)
        for jj in range(2):  # We iterate over 2 because each i constraint has an upper and lower bound (i.e 2 constraints)
            for ii in range(self.n):
                if jj == 0: # Upper bound constraint
                    zcbf_p[ii] = self.beta.eval(self.b_bar(ii,q[ii],v[ii])) 
                    Lf_b2[ii] = -self.gamma*self.alpha.diff(self.h_bar(ii,q[ii]))*v[ii]
                if jj == 1: # Lower bound constraint
                    zcbf_p[ii+self.n] = self.beta.eval(self.b_ubar(ii,q[ii],v[ii])) 
                    Lf_b2[ii+self.n] = self.gamma*self.alpha.diff(self.h_ubar(ii,q[ii]))*v[ii]

        ones_2n = np.ones(2*self.n)
        zcbf_b = -Lf_b1 - Lf_b2 - self.nu*zcbf_p + ones_2n*self.eta

        # Define QP parameters for qp solver (CVXOPT)
        qp_P = cvxopt.matrix(np.eye(self.n))
        qp_q = cvxopt.matrix(-u_nom)
        qp_C_zcbf = -Lg_b
        qp_C_umax = np.eye(self.n)
        qp_C_umin = -np.eye(self.n)
        qp_G = cvxopt.matrix(np.concatenate((qp_C_zcbf, qp_C_umax, qp_C_umin)))
        qp_h = cvxopt.matrix(np.concatenate((-zcbf_b, np.array(self.u_max), -np.array(self.u_min))))

        cvxopt.solvers.options['show_progress'] = False # mute optimization output
        cvxopt.solvers.options['maxiters'] = 500 # increase max iteration number
        cvxopt.solvers.options['abstol'] = 1e-4
        cvxopt.solvers.options['reltol'] = 1e-5
        try:
            solv_star = cvxopt.solvers.qp(qp_P,qp_q,qp_G,qp_h)
        except:
            print('cvxopt solver failed:')
            print('time: ' + str(t))
            print('(q,v) in region ' + str(self.determine_region_H(q,v) ))
            print(' ubar = ' + str(self.compute_ubar_control(q,v)))
            print('qp_P: ' + str(qp_P))
            print('qp_q: ' + str(qp_q))
            print('qp_G: ' + str(qp_G))
            print('qp_h: ' + str(qp_h))
            raise TypeError('no control input found')
            #cvxopt.solvers.options['show_progress'] = True
            #solv_star = cvxopt.solvers.qp(qp_P,qp_q,qp_G,qp_h)

        u_star = np.array(solv_star['x'])# extract decision variables of optimization problem

        return u_star.flatten()  

    def closed_loop_dynamics(self,x,t,k,uk=None):

        # Extract states
        q = x[:self.n]
        v = x[self.n:]

        # Compute system dynamics
        try:
            f_x, g_x = self.sys_dyn.compute_system_dynamics(q,v,t)
        except:
            raise TypeError('No system dynamics defined!')

        # Compute control 
        if k == 'ZCBF_control':
            u = self.compute_zcbf_control(q,v,f_x,g_x,t)
        elif k == 'nom_control':
            u = self.sys_dyn.nominal_control(q,v,t)
        elif k == 'ZCBF_control_sampled':
            u = uk
        else:
            u = np.zeros(self.n)

        # Compute state derivatives
        qdot = v
        vdot = f_x + g_x.dot(u)
        
        xdot = np.concatenate((qdot,vdot))

        # Store control input
        #self.store_inputs(t,u)
        return xdot

    def compute_zcbf_u_trajectory(self,t,q_traj,v_traj,k):
        '''Recompute control values used in odeint'''

        u_traj = np.zeros((len(t),self.n))

        for ii in range(len(t)):
            ti = t[ii]
            qi = q_traj[ii,:]
            vi = v_traj[ii,:]
            if k == 'ZCBF_control':
                f_x, g_x = self.sys_dyn.compute_system_dynamics(qi,vi,ti)
                u_traj[ii,:]= self.compute_zcbf_control(qi,vi,f_x,g_x,ti)
            if k == 'nom_control':
                u_traj[ii,:]= self.sys_dyn.nominal_control(qi,vi,ti)

        return u_traj


    def run_simulation(self,t,x0,cont_type,rtol=None,atol=None,hmax=0.0,cont=True,T=None):
        '''Run odeint for ZCBF simulations. cont used to define if control inut is continuous
        or discrete (i.e. zero-order hold) with T sampling time'''

        print('------------')
        print('Running simulation...')

        # Run odeint depending on if control is continuous or discrete
        # For continuous control, run odeint as usual and then compute control from resulting state trajectory
        if cont:
            print('(cont time)')
            y = odeint(self.closed_loop_dynamics,x0,t,args=(cont_type,),rtol=rtol,atol=atol,hmax=hmax)
            u = self.compute_zcbf_u_trajectory(t,y[:,0:2],y[:,2:4],cont_type)
            tout = t

        else: # For discrete control, run odeint between sampling times. Control is computed at each new sampling time and held constant between
            print('(sampled-data)')
            k = 0 # initialize sampling time counter
            x0_k = x0 # Initialize initial condition for k = 0
            y = np.zeros((1,2*self.n)) # initialize state trajectory
            y[0,:] = x0
            N_c = (t[-1] - t[0])/(t[1] - t[0]) + 1.0 # Compute number of time data points specified for the integrator
            N_s = (t[-1] - t[0])/T # COmpute number of time data points for sampled time
            Nk = int(round(N_c/N_s)) # Number of data points for time tk = [kT,...(k+1)T]
            if Nk < 2:
                Nk = 2
            #u = np.zeros((1,self.n))    # Setup control list
            tout = [t[0]]
            while k*T < t[-1]:

                # Extract sampled states and compute control
                qk = x0_k[0:2]
                vk = x0_k[2:4]
                f_x, g_x = self.sys_dyn.compute_system_dynamics(qk,vk,k*T)
                uki = self.compute_zcbf_control(qk,vk,f_x,g_x,k*T)

                # Compute time array and repeat control values over time array
                tk = np.linspace(k*T, (k+1)*T, Nk)
                uk = np.tile(uki,(Nk,1))
                if k == 0:
                    u = uk
                else:
                    u = np.concatenate((u,uk[1:,:]))

                # Run odeint on constant control value uki
                yk = odeint(self.closed_loop_dynamics,x0_k,tk,args=(cont_type,uki),rtol=rtol,atol=atol,hmax=hmax)
                y = np.concatenate((y,yk[1:,:]))
                tout.extend(tk[1:].tolist())

                # Increment discrete counter and reset initial condition
                k += 1
                x0_k = yk[-1,:]

                # print(uki)
                # print(tk)
                # print(uk)
                # if k > 1:
                #     print(tout)
                #     #print(y)
                #     print(u)
                #     raise(SystemError)

        return tout,y,u

    def compute_zcbf_parameters(self,sim_control,dq = 0.5,T_sample=0.0,eps_frac = 0.5,nu_frac=0.5):
        '''Compute safe-by-design zcbf parameters'''

        # Compute epsilon and gamma2_star, which quantifies how much excess control authority the system has
        # and the max velocity allowed in the safe set respectively so that the control can enforce foward invariance
        # To compute epsilon and gamm2_star, we do a grid search over Qdelta for all i 

        qmin_delta = self.q_min - self.delta*np.ones(self.n)
        qmax_delta = self.q_max + self.delta*np.ones(self.n)
        dq_list = [dq for ii in range(self.n)]
        
        # Compute epsilon parameter
        max_epsilon = two_dof_module.grid_loop(qmin_delta.tolist(), qmax_delta.tolist(), dq_list, self.compute_eps_function, min)
        self.epsilon = eps_frac*max_epsilon
        print('-----------------')
        print('Epsilon value computed: as ' + str(eps_frac)+'*max epsilon')
        print('max epsilon: ' + str(max_epsilon))
        print('epsilon: ' + str(self.epsilon) )

        # Compute gamma2_star parameter
        #self.gamma2_star = self.compute_grid_search(self.compute_gamma2_star_function,dq)
        #self.gamma2_star = two_dof_module.compute_grid_search(self.compute_gamma2_star_function,dq,self.q_min - self.delta*np.ones(self.n),self.q_max + self.delta*np.ones(self.n))
        self.gamma2_star = two_dof_module.grid_loop(qmin_delta.tolist(), qmax_delta.tolist(), dq_list, self.compute_gamma2_star_function, min)
        print('-----------------')
        print('gamma2_star value computed:')
        print('gamma2_star = ' + str(self.gamma2_star))

        # Compute gamma1 star parameter
        a = self.alpha.eval(2*self.delta + max(np.absolute(np.array(self.q_max) - np.array(self.q_min))))
        self.gamma1_star = min(self.v_max)/a
        print('-----------------')
        print('gamma1_star value computed:')
        print('gamma1_star = ' + str(self.gamma1_star))

        # Compute gamma3 star parameter
        L = self.alpha.eval_L()
        self.gamma3_star = math.sqrt(self.epsilon/(L*a))
        print('-----------------')
        print('gamma3_star value computed:')
        print('gamma3_star = ' + str(self.gamma3_star))

        # Choose gamma
        self.gamma = min(self.gamma, self.gamma1_star, self.gamma2_star, self.gamma3_star)
        print('-----------------')
        print('gamma value chosen as min(gamma(original), gamma1_star, gamma2_star, gamma3_star)')
        print('gamma = ' + str(self.gamma))

        # Compute delta_star
        self.delta_star = self.compute_delta_star()
        self.delta = min(self.delta, self.delta_star)

        print('-----------------')
        print('delta_star value computed')
        print('delta_star = ' + str(self.delta_star))
        print('delta = min(delta,delta_star):')
        print('delta = ' + str(self.delta))

        # Compute zeta_delta, rho_delta and v_bar
        zeta_delta = self.alpha.eval_zeta_delta(self.delta,self.gamma)
        rho_delta = self.alpha.eval_rho_delta(self.delta,self.gamma,self.q_max,self.q_min,self.rho)
        self.v_bar = self.gamma*a

        # Compute nu1_star, nu2_star, and choose nu
        self.nu1_star = ( self.gamma**2*L*a )/( self.beta.eval(rho_delta) )
        self.nu2_star = self.epsilon/( abs(self.beta.eval(zeta_delta) ))
        

        print('-----------------')
        print('nu1_star,nu2_star values computed:')
        print('nu1_star = ' + str(self.nu1_star))
        print('nu2_star = ' + str(self.nu2_star))
        if (self.nu > self.nu2_star) or (self.nu <= self.nu1_star):
            if self.delta > 0:
                print('Original nu inadequate, choosing nu = ' +str(nu_frac)+'*(nu2_star -nu1_star) + nu1_star for delta > 0')
                self.nu = nu_frac*self.nu2_star + self.nu1_star*(1.0 - nu_frac)
            else:
                print('Original nu indadequate, choosing nu = (1 - '+ str(nu_frac)+') nu1_star for delta = 0: ')
                self.nu = (1.0 - nu_frac)*self.nu1_star
        else:
            print('Original nu in (nu1_star, nu2star], so no change:')
        print('nu = ' + str(self.nu))

        # Compute eta_star, eta
        self.eta_star = 0.5*( self.nu*self.beta.eval(rho_delta) - self.gamma**2*L*a  )
        print('-----------------')
        print('eta_star value computed:')
        print('eta_star = ' + str(self.eta_star))
        if sim_control == 'discrete_time':
            
            print('eta_bar = ' + str(self.eta))
            if self.eta_star < 0:
                print('eta_star = ' + str(self.eta_star))
                raise TypeError('In computing eta_star, eta_star <0 was encountered!')

            if self.eta > self.eta_star:
                print('Encountered eta_bar > eta_star, setting eta_bar = eta_star')
                self.eta = self.eta_star
                print('eta_bar = ' + str(self.eta))
                
        

            #Compute c2, lipschitz constant of beta. Requires largest value of b depending on type of function
            bmax = 0.0
            for ii in range(self.n):
                bmax = max(bmax, self.b_bar(ii, self.q_min[ii] - self.delta, -self.v_max[ii] ) )
            c2 = self.beta.eval_L(bmax)

            # Compute c4 which is simply the max possible control torque i.e. ||u_max||_infty
            c4 = np.linalg.norm(self.u_max, ord=np.inf)

            # Compute c5 which is the upper bound on the dynamics
            #kf = max(np.diag(self.sys_dyn.F))
            #c5 = self.sys_dyn.km2*( self.sys_dyn.kc*self.v_bar**2 + kf*self.v_bar + self.sys_dyn.kg + c4 )

            # Compute eta(T):
            eta_num = ( self.sys_dyn.c1 + c2 + self.sys_dyn.c3*c4 )*self.sys_dyn.c5
            eta_den = self.sys_dyn.c1 + c2*c4 
            eta_T = eta_num*math.expm1( eta_den*T_sample )/eta_den
            print('eta(T_sample): '+ str(eta_T))

            #eta_lin_T = eta_num*T_sample
            #print('eta_lin(T_sample): ' + str(eta_lin_T))
            print('requirement: eta(T) <= eta_bar <= eta_star')

            
            print('T_sample = ' + str(math.log1p(eta_T*eta_den/eta_num)/eta_den ))

            if eta_T > self.eta:
                eta_inv = math.log1p(self.eta*eta_den/eta_num)/eta_den
                raise TypeError('Sampling time not fast enough! eta(T_sample) > eta_bar. Choose T_sample <= ' + str(eta_inv))

            print('-----------------')
            print('eta value chosen: eta = eta(T)')
            self.eta = eta_T
            print('eta = ' + str(self.eta))

        else: # if cont_time, no need for eta, set to zero
            self.eta = 0.0
            print('Since simulation in cont. time: eta set to zero, eta = ' + str(self.eta))

    def compute_delta_star(self, N_delta = 100):
        '''Compute delta star parameter'''

        delta_range = np.linspace(0.0,self.delta,N_delta)
        delta_stop = False

        # Iterate over delta_range checking until the following condition fails: |beta(zeta_delta)| <= beta(rho_delta)
        # delta_star is the maximum value for which the condition holds
        for delta_iter in delta_range:
            zeta_delta = self.alpha.eval_zeta_delta(delta_iter,self.gamma)
            rho_delta = self.alpha.eval_rho_delta(delta_iter,self.gamma,self.q_max,self.q_min,self.rho)


            if ( abs(self.beta.eval( zeta_delta )) <= self.beta.eval( rho_delta ) ) and (delta_stop != True):
                delta_star = delta_iter
            else:
                delta_stop = True



        return delta_star
        

    def compute_eps_function(self,ii,q,dq,data):
        '''Compute function to minimize in order to compute epsilon'''
        g = self.sys_dyn.g_eval(q)
        M = self.sys_dyn.M_eval(q)
        ei = np.zeros((1,self.n))
        ei[0,ii] = 1.0


        # Function to minimize to compute epsilon
        fi = (self.u_max[ii] - abs(g[ii]))/( max(np.absolute(np.matmul(ei,M).flatten())) ) - self.eta

        if fi <= 0: # check if infeasible delta/eta/umax/umin were chosen
            print('epsilon <= ')
            print(fi)
            raise TypeError('In computing epsilon, epsilon <=0 was encountered!. Choose eta_bar smaller or increase input constraints')

        return fi

    def compute_gamma2_star_function(self,ii,q, dq, data):
        '''Compute function to minimize in order to compute gamma2_star'''
        g = self.sys_dyn.g_eval(q)
        M = self.sys_dyn.M_eval(q)
        Mi_max = max(np.absolute(M[ii,:]))
        a = self.alpha.eval(2.0*self.delta + max(np.absolute(np.array(self.q_max) - np.array(self.q_min))))
        y = self.y_eval(q)

        di = self.sys_dyn.F[ii,ii]/( Mi_max*y + self.sys_dyn.kc*a )
        ci = (abs(g[ii]) + (self.epsilon + self.eta)*Mi_max- self.u_max[ii] )/( Mi_max*y*a + self.sys_dyn.kc*a**2 )

        gi = 0.5*( -di + math.sqrt(di**2 - 4.0*ci) )

        if gi <= 0:
            print('gamma2_star <=')
            print(gi)
            print(Mi_max)
            print(y)
            print(a)
            raise TypeError('In computing gamma2_star, gamma2_star <=0 was encountered!')

        return gi

    def y_eval(self,q):
        '''compute function y = max dalpha/dhbar, dalphadhubar over all i = 1,..,n'''
        y = 0.0
        for ii in range(self.n):
            y = max(self.alpha.diff(self.h_bar(ii,q[ii]) ), self.alpha.diff(self.h_ubar(ii,q[ii])))

        return y

    



