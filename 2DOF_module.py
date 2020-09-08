import numpy as np
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt


class 2DOFSystem():
    def __init__(self,m_list,l_list,f_list,I_list,grav):
        '''Initialize class defining 2DOF system dynamics'''

        # Store model parameters
        self.n_links = 2
        if (self.n_links == len(m_list) and self.n_links == len(l_list) 
                and self.n_links== len(f_list) and self.n_links == len(I_list)):
            self.m_list = m_list
            self.l_list = l_list
            self.f_list = f_list
            self.I_list = I_list
        else:
            raise TypeError('m_list...,I_list not of same length. (2DOFSystem class)') 

        self.grav = grav
        
    def M_eval(self,q):
        '''Compute inertia matrix'''

        # Extract parameters
        m0 = self.m_list[0]
        m1 = self.m_list[1]
        l0 = self.l_list[0]
        l1 = self.l_list[1]
        I0 = self.I_list[0]
        I1 = self.I_list[1]

        M = np.zeros(n_links)
        M[0,0] = m0*(l0/2)**2 + m1*(l0**2 + (l1/2)**2 + 2*l0*(l1/2)**2 + 2*l0*(l1/2)*cos(q[1])) + I0 + I1;
        M[0,1] = m1*( (l1/2)**2 + l0*(l1/2)*cos(q[1]) ) + I1;
        M[1,0] = M[0,1];
        M[1,1] = m1*(l1/2)**2 + I1;

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

        C = np.zeros(n_links)
        C[0,0] = -m1*l0*(l1/2)*sin(q[1]);
        C[0,1] = -m1*l0*(l1/2)*sin(q[1]);
        C[1,0] =  m1*l0*(l1/2)*sin(q[1]);
        C[1,1] = 0;

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

        g = np.zeros(n_links,1)
        g[0] = (m0*(l0/2)+m1*l0)*self.grav*sin(q[0]) + m1*self.grav*(l1/2)*sin(q[0]+q[1]);
        g[1] =  m1*self.grav*(l1/2)*sin(q[0]+q[1]);

        return g