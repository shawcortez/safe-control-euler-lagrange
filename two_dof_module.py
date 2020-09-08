import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt



class TwoDOFSystem():
    def __init__(self,m_list,l_list,f_list,I_list,grav,gains,ref, nom_control_type, kc=None, kg=None, km2=None, c3 = None, c1 = None, c5 = None):
        '''Initialize class defining 2DOF system dynamics'''

        # Store model parameters
        self.n_joints = 2
        if (self.n_joints == len(m_list) and self.n_joints == len(l_list) 
                and self.n_joints== len(f_list) and self.n_joints == len(I_list)):
            self.m_list = m_list
            self.l_list = l_list
            self.f_list = f_list
            self.I_list = I_list
        else:
            raise TypeError('m_list...,I_list not of same length. (2DOFSystem class)') 

        self.grav = grav
        self.F = np.diag(np.array(f_list))
        self.gains = gains
        self.ref = ref
        self.nom_control_type = nom_control_type
        self.kc = kc
        self.kg = kg
        self.km2 = km2
        self.c3 = c3
        self.c1 = c1
        self.c5 = c5

    def M_eval(self,q):
        '''Compute inertia matrix'''

        # Extract parameters
        m0 = self.m_list[0]
        m1 = self.m_list[1]
        l0 = self.l_list[0]
        l1 = self.l_list[1]
        I0 = self.I_list[0]
        I1 = self.I_list[1]

        M = np.zeros((self.n_joints,self.n_joints))

        M[0,0] = m0*(l0/2.0)**2 + m1*(l0**2 + (l1/2.0)**2 + 2.0*l0*(l1/2.0)**2 + 2.0*l0*(l1/2.0)*np.cos(q[1])) + I0 + I1;
        M[0,1] = m1*( (l1/2.0)**2 + l0*(l1/2.0)*np.cos(q[1]) ) + I1;
        M[1,0] = M[0,1];
        M[1,1] = m1*(l1/2.0)**2 + I1;

        return M

    def C_eval(self,q,v):
        '''Compute Corilois/centrifugal matrix'''

        # Extract parameters
        m0 = self.m_list[0]
        m1 = self.m_list[1]
        l0 = self.l_list[0]
        l1 = self.l_list[1]
        I0 = self.I_list[0]
        I1 = self.I_list[1]

        C = np.zeros((self.n_joints,self.n_joints))
        C[0,0] = -m1*l0*(l1/2.0)*np.sin(q[1])*v[1];
        C[0,1] = -m1*l0*(l1/2.0)*np.sin(q[1])*(v[0] + v[1]);
        C[1,0] =  m1*l0*(l1/2.0)*np.sin(q[1])*v[0];
        C[1,1] = 0.0;

        return C

    def g_eval(self,q):
        '''Compute gravity torque'''

        # Extract parameters
        m0 = self.m_list[0]
        m1 = self.m_list[1]
        l0 = self.l_list[0]
        l1 = self.l_list[1]
        I0 = self.I_list[0]
        I1 = self.I_list[1]

        g = np.zeros(self.n_joints)
        g[0] = (m0*(l0/2.0)+m1*l0)*self.grav*np.sin(q[0]) + m1*self.grav*(l1/2.0)*np.sin(q[0]+q[1]);
        g[1] =  m1*self.grav*(l1/2.0)*np.sin(q[0]+q[1]);

        return g

    def nominal_control(self,q,v,t):
        '''Compute a pre-defined nominal control law for the system'''
        if self.nom_control_type == 'PD':
            u_nom = -self.gains['Kp']*q - self.gains['Kv']*v

        if self.nom_control_type == 'ComputedTorque':
            u_nom = self.computed_torque_control(q,v,t)


        return u_nom

    def computed_torque_control(self,q,v,t):
        '''Computed torque control law from Spong for nominal controller'''

        # Compute reference and error states
        r, dr, ddr = self.computed_torque_reference(t)
        e = q - r
        edot = v - dr

        # Compute dynamic terms
        M = self.M_eval(q)
        C = self.C_eval(q,v)
        g = self.g_eval(q)

        return M.dot(ddr - self.gains['Kv']*edot - self.gains['Kp']*e) + C.dot(v) + g



    def computed_torque_reference(self,t):
        '''Compute reference trajectory for nominal controller'''

        r0 = self.ref['r0_amp']*np.sin(self.ref['r0_omega']*t)
        r0_dot = self.ref['r0_amp']*self.ref['r0_omega']*np.cos(self.ref['r0_omega']*t)
        r0_ddot = -self.ref['r0_amp']*self.ref['r0_omega']**2*np.sin(self.ref['r0_omega']*t)

        r1 = self.ref['r1_amp']*np.sin(self.ref['r1_omega']*t) + self.ref['r1_b']
        r1_dot = self.ref['r1_amp']*self.ref['r1_omega']*np.cos(self.ref['r1_omega']*t)
        r1_ddot = -self.ref['r1_amp']*self.ref['r1_omega']**2*np.sin(self.ref['r1_omega']*t)

        ref = np.array([r0, r1])
        ref_dot = np.array([r0_dot, r1_dot])
        ref_ddot = np.array([r0_ddot, r1_ddot])
        return ref, ref_dot, ref_ddot

    def compute_system_dynamics(self,q,v,t):
        '''Compute the system dynamics f_x, g_x '''

        # Compute system dynamics
        M = self.M_eval(q)
        C = self.C_eval(q,v)
        g = self.g_eval(q)
        M_inv = np.linalg.inv(M)
        f_x = M_inv.dot( -C.dot(v) - self.F.dot(v) + g)
        g_x = M_inv

        return f_x, g_x;

    def compute_system_bounds(self,q_min,q_max,v_min, v_max,dq_ = 0.1, dv_ = 0.1):
        '''Compute bounds on the system dynamics, km2, kg, kc, kf'''

        # Construct delta lists for q and v if needed
        if isinstance(dq_,list):
            dq = dq_
        else:
            dq = [dq_ for q in q_min]

        if isinstance(dv_,list):
            dv = dv_
        else:
            dv = [dv_ for v in v_min]

        # Compute km2, kc, kg
        km2 = grid_loop(q_min, q_max, dq, self.compute_km2_func, max)
        print('km2 computed: ' + str(km2))

        #raise(SystemExit)

        kc = grid_loop(q_min+v_min, q_max+v_max, dq+dv, self.compute_kc_func,max)
        print('kc computed: ' + str(kc))

        #raise(SystemExit)

        kg = grid_loop(q_min, q_max, dq, self.compute_kg_func, max)
        print('kg computed: ' + str(kg))

        #raise(SystemExit)

        # Compute Lipschitz bounds on inv(M)
        c3 = grid_loop(q_min, q_max, dq, self.compute_L_M_func, max)
        print('c3 computed: ' + str(c3))

        # Compute system bounds i.e. || inv(M)*( - Cv - Fv - g ) ||_infty + ||inv(M)||_infty ||u||_infty
        c5 = grid_loop(q_min+v_min, q_max+v_max, dq+dv, self.compute_max_dynamics_func, max)
        print('c5 computed: ' + str(c5))

        # Compute Lipschitz bounds on inv(M)*( - Cv - Fv - g )
        c1 = grid_loop(q_min+v_min, q_max+v_max, dq+dv, self.compute_L_dynamics_func, max)
        print('c1 computed: ' + str(c1))


        return

    def compute_km2_func(self,ii,q,dq,data):
        '''function to minimize to compute km2. Equates to finding the max eigenvalue of inv(M)'''
        M = self.M_eval(q)

        #print(np.linalg.eigvals(np.linalg.inv(M)))
        #print(q)

        return max(np.linalg.eigvals(np.linalg.inv(M)))

    def compute_kc_func(self,ii,x,dx,data):
        '''function to minimize to compute kc. Equates to finding max value of || C(q,1) || (as an approximation)'''
        q = x[0:2]
        v = x[2:4]
        
        # Compute system dynamics
        M = self.M_eval(q)
        C = self.C_eval(q,v)
        g = self.g_eval(q)

        if np.linalg.norm(v,ord=np.inf) != 0.0:
            kc_i = np.linalg.norm(C,ord=np.inf)/ np.linalg.norm(v,ord=np.inf)
        else:
            kc_i = 0.0

        #print(kc_i)
        #print(C)
        #print(np.linalg.norm(C,ord=np.inf))
        #print(np.linalg.norm(C))

        return kc_i
        #return -np.linalg.norm(C)

    def compute_kg_func(self,ii,q, dq,data):
        '''function to minimize to compute kg. Equates to finding max value of || g(q) ||'''
        g = self.g_eval(q)

        return np.linalg.norm(g,ord=np.inf)

    def compute_L_dynamics_func(self,ii,x,dx ,data):
        '''computes Lipschitz constant of dynamics by computing the slopt of the dynamics for all grid points 
        around the current x'''

        q = x[0:2]
        v = x[2:4]
        xmin_array = np.array(x) - np.array(dx)
        xmax_array = np.array(x) + np.array(dx)
        xmin = xmin_array.tolist()
        xmax = xmax_array.tolist()

        return grid_loop(xmin,xmax,dx,self.compute_L_dynamics_Lipschitz,max, data = x)

    def compute_L_dynamics_Lipschitz(self,ii,x,dx,data):
        '''computes the slop of the dynamics over dq for the current x w.r.t the centerpoint, data'''
        q_eps = data[0:2]
        v_eps = data[2:4]
        qi = x[0:2]
        vi = x[2:4]

        M_inv = np.linalg.inv(self.M_eval(q_eps))
        C = self.C_eval(q_eps,v_eps)
        g = self.g_eval(q_eps)
        dyn_func = np.matmul(M_inv, -np.matmul(C,v_eps) - np.matmul(self.F, v_eps) - g)

        Mi_inv = np.linalg.inv(self.M_eval(qi))
        Ci = self.C_eval(qi,vi)
        gi = self.g_eval(qi)
        dyn_func_i = np.matmul(Mi_inv, -np.matmul(Ci,vi) - np.matmul(self.F, vi) - gi)

        if np.linalg.norm(np.array(data) - np.array(x)) == 0.0:
            fi = 0.0
        else:
            fi = (np.linalg.norm(dyn_func_i - dyn_func) )/ np.linalg.norm(np.array(data) - np.array(x))

        return fi


    def compute_L_M_func(self,ii,q,dq ,data):
        '''Computes Lipschitz constant of inv(M). Uses current q from outer loop to check all neighboring
        grid points to then compute the slope of inv(M) w.r.t all neighboring grid points of data'''

        qmin_array = np.array(q) - np.array(dq)
        qmax_array = np.array(q) + np.array(dq)
        qmin = qmin_array.tolist()
        qmax = qmax_array.tolist()
        
        return grid_loop(qmin,qmax,dq,self.compute_L_M_Lipschitz,max, data = q)

        

    def compute_L_M_Lipschitz(self,ii,q,dq,data):
        '''computes slope: (inv(M) - inv(M))/dq for Lipschitz computation'''

        M_inv = np.linalg.inv(self.M_eval(data))

        Mi_inv = np.linalg.inv(self.M_eval(q))

        if np.linalg.norm(np.array(data) - np.array(q)) == 0.0:
            fi = 0.0
        else:
            fi = (np.linalg.norm(Mi_inv - M_inv) )/ np.linalg.norm(np.array(data) - np.array(q))

        return fi

    def compute_max_dynamics_func(self,ii,x,dx ,data):
        '''Computes max of inv(M)(-Cv - Fv - g + u) . Uses current q from outer loop to check all neighboring
        grid points to then compute the slope of inv(M) w.r.t all neighboring grid points of data'''

        q = np.array(x[0:2])
        v = np.array(x[2:4])

        M_inv = np.linalg.inv(self.M_eval(q))
        C = self.C_eval(q,v)
        g= self.g_eval(q)

        fi = np.linalg.norm(np.matmul(M_inv, -np.matmul(C,v) - np.matmul(self.F, v)  -g), ord=np.inf) + np.linalg.norm(M_inv, ord=np.inf)*np.linalg.norm(np.array(self.u_max),ord=np.inf)
        
        return fi


def grid_loop(xmin,xmax,dx,func_eval,func_loop=min, data = None):
    '''Peform grid search over x variable. Start at xmin and increment each component of x by dx resp.
    At each loop, we evaluate some function, func_eval, and evaluate a function, func_loop, over each 
    evaluateion of func_eval in the loop.'''

    n_loops = len(xmin)

    return loop_me(n_loops-1,xmin[:],xmax,dx,func_eval,func_loop,data = data)


def loop_me(ii,x,xmax,dx,func_eval,func_loop=min, fi=None, data = None):
    '''Perform recursive grid loop over each element of x and evaluate func_eval for each x. Then evaluate
    func_loop over each func_eval. Return the value of func_loop evaluated over all loops.'''

    #print('x starting loop:' + str(x))
    #print('ii starting loop:' + str(ii))

    while x[ii] <= xmax[ii]:

        if ii > 0:
            
            #print('enter loop')
            fi = loop_me(ii-1,x[:],xmax,dx,func_eval,func_loop, fi, data)

            
            #print('x after loop: ' + str(x))

        else:
            #print(x)
            
            if fi == None:
                #print('fi = None')
                fi = func_eval(ii,x,dx, data)
            else:
                fi = func_loop(fi, func_eval(ii,x,dx,data))

        # increment step
        x[ii]  += dx[ii]

        
        #print(fi)
        #print(ii)

    #print('exit loop')
    #print(x)

    return fi

def test_func(ii,x,dx):
    return np.linalg.norm(np.array(x))

#     def compute_grid_search(func_min,dq,qmin,qmax):
#     '''Compute grid search over all Qdelta and i = 1,2 to minimize the given func_min'''

#     print('Running grid search')

#     q_eps = 1.0*np.array(qmin) # initialize grid search
#     fmin = None # initialize function to minimize
#     n = len(qmin)


#     kk = 0 # outerloop counter
#     while q_eps[1] <= qmax[1]: # Iterate over second dimension

#         # Reset middleloop counter on first dimension
#         q_eps[0] = qmin[0]
#         jj = 0

#         while q_eps[0] <= qmax[0]: # Iterate over first dimension
            
#             for ii in range(n): # Compute fmin,and gmin over every ii on current q_eps

#                 # Function to minimize to compute epsilon
#                 fi = func_min(ii,q_eps)

#                 # Store values of fi, searching for minimum value over all q, then all ii
#                 if fmin == None: # initialize fmin value
#                     fmin = fi  
#                 else:            # ensure fmin is the min of all fi values (over ii and q_eps)
#                     fmin = min(fi,fmin)


#                 # if kk == 0:
#                 #     print('Ending sim')
#                 #     raise(SystemExit)

#             # increment middle loop counter
#             q_eps[0] += dq
#             jj += 1

#         # increment outerloop counter
#         q_eps[1] += dq
#         kk += 1

#     return fmin

# def compute_full_grid_search(func_min,dq,dv,qmin,qmax,vmin,vmax):
#     '''Compute grid search over all Qdelta and i = 1,2 to minimize the given func_min'''

#     print('Running full grid search')

#     x_eps = np.concatenate((qmin, vmin),axis=None) # initialize grid search
#     fmin = None # initialize function to minimize
#     n = len(x_eps)


#     kk = 0 # outerloop counter
#     while x_eps[3] <= vmax[1]: # Iterate over second dimension

#         # Reset middleloop counter on first dimension
#         x_eps[2] = vmin[0]
#         jj = 0

#         while x_eps[2] <= vmax[0]: # Iterate over first dimension

#             # Reset
#             x_eps[1] = qmin[1]

#             while x_eps[1] <= qmax[1]:

#                 # Reset
#                 x_eps[0] = qmin[0]

#                 while x_eps[0] <= qmax[0]:

#                     #print('x_eps: ' + str(x_eps))
            
#                     for ii in range(n): # Compute fmin,and gmin over every ii on current q_eps

#                         # Function to minimize to compute epsilon
#                         fi = func_min(ii,x_eps)

#                         # Store values of fi, searching for minimum value over all q, then all ii
#                         if fmin == None: # initialize fmin value
#                             fmin = fi  
#                         else:            # ensure fmin is the min of all fi values (over ii and q_eps)
#                             fmin = min(fi,fmin)


#                 # if kk == 0:
#                 #     print('Ending sim')
#                 #     raise(SystemExit)

                
#                     # increment middle loop counter
#                     x_eps[0] += dq

#                 # increment middle loop counter
#                 x_eps[1] += dq

#             # increment middle loop counter
#             x_eps[2] += dq

#         # increment outerloop counter
#         x_eps[3] += dq
#         kk += 1

#     return fmin
