import numpy as np
import math
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import two_dof_module
import cvxopt
from scipy.integrate import odeint

class NonlinearConstraints():
    def __init__(self, type, params):

        # Define constraint type
        self.type = constraint_type

        # Store constraint parameters
        self.a = params['a'] 
        P_array = np.array(params['P'])
        n = len(P_array)
        self.P = np.reshape(P_array, (int(math.sqrt(n)),int(math.sqrt(n))) )
        self.qr = np.array(params['qr'])

    def eval(self,x):
        '''Evalute the constraint'''

        if self.type == 'ellipsoid':
            return self.a - 0.5*np.dot( x - self.qr, self.P.dot(x-self.qr))

        if self.type == 'planar':
            return self.qr.T.dot(x) + self.a

    def eval_gradient(self, x):
        '''Evaluate the constraint gradient at x'''

        if self.type == 'ellipsoid':
            return -0.5*(self.P + self.P.T).dot( q-self.qr )

        if self.type == 'planar':
            return self.qr

    def eval_hessian(self, x):
        '''Evaluate the constraint hessian at x'''

        if self.type == 'ellipsoid':
            return -0.5*(self.P + self.P.T)

        if self.type == 'planar':
            return np.zeros((len(x), len(x)))


class NonlinearTransformation():
    def __init__(self, dynamics, constraint_type_list, params_list):

        # Store the unstransformed dynamical system
        self.dynamics = dynamics
        self.n = dynamics.n_joints # number of position states in the system

        # Check that the number of constraints is equal to the number of position states
        if (len(constraint_type_list) != self.n) or (len(params_list != self.n)):
            raise TypeError('Not enough constraints defined. Need n constraints') 

        # Construct constraints
        for ii, ci_type in enumerate(constraint_type_list):
            self.c[ii] = NonlinearConstraints(ci_type,params_list[ii])

    def eval(self, x):
        '''Transform the given state of the original dynamics'''

        return np.array([self.c[ii].eval(x) for ii in range(self.n)]  )

    def eval_gradient(self, x):
        '''Compute the gradient of the constraint transformation'''

        return np.array([self.c[ii].eval_gradient(x) for ii in range(self.n)]  )

    def eval_hessian(self, x):
        '''Compute the hessian of the constraint transformation'''

        # TO DO: write in tensors 







        