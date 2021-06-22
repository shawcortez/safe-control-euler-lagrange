import numpy as np
import math
import matplotlib
matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
import two_dof_module
import cvxopt
from scipy.integrate import odeint

class NonlinearConstraints():
    def __init__(self, constraint_type, params):

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
            return -0.5*(self.P + self.P.T).dot( x-self.qr )

        if self.type == 'planar':
            return self.qr

    def eval_hessian(self, x):
        '''Evaluate the constraint hessian at x'''

        if self.type == 'ellipsoid':
            return -0.5*(self.P + self.P.T)
            #return -0.0*(self.P + self.P.T) #CHANGE THIS IS FOR DEBUGGING

        if self.type == 'planar':
            return np.zeros((len(x), len(x)))

class IDConstraints():
    def __init__(self,ii,n):
        '''Define identity function for with state of size n'''
        self.ii = ii
        self.n = n
        self.ei = np.eye(n)[:,ii]

    def eval(self, x):
    
        return x[self.ii]

    def eval_gradient(self, x):

        return self.ei

    def eval_hessian(self, x):

        return np.zeros((self.n, self.n))


class NonlinearTransformation():
    def __init__(self, n, constraint_type_list=None, params_list=None):

        # Store the unstransformed dynamical system
        self.n = n # number of position states in the system

        try: # If constraint types non-empty
            # Check that the number of constraints is equal to the number of position states
            if (len(constraint_type_list) != self.n) or (len(params_list) != self.n):
                raise TypeError('Not enough constraints defined. Need n constraints') 

            # Construct constraints
            self.c = []
            for ii, ci_type in enumerate(constraint_type_list):
                self.c.append( NonlinearConstraints(ci_type,params_list[ii]) )


            self.id_bool = False

        except: 
            # If no constraints defined, use identity functions
            # Construct constraints
            self.c = []
            for ii in range(self.n):
                self.c.append( IDConstraints(ii,self.n)  )

            self.id_bool = True

       

    def eval(self, x):
        '''Transform the given state of the original dynamics'''

        return np.array([self.c[ii].eval(x) for ii in range(self.n)]  )

    def eval_gradient(self, x):
        '''Compute the gradient of the constraint transformation'''

        # Q: Should the gradient term be transposed??
        return np.array([self.c[ii].eval_gradient(x).T for ii in range(self.n)]  )

    def eval_hessian(self, x):
        '''Compute the hessian of the constraint transformation'''

        # Q: Should the hessian term be transposed as well??
        return np.array([self.c[ii].eval_hessian(x) for ii in range(self.n)]  )

    def compute_hessian_infnorm(self, x):
        '''Compute infinitey norm of hessian'''

        norm_i = [np.linalg.norm(self.c[ii].eval_hessian(x),ord=np.inf) for ii in range(self.n)]
        return max(norm_i)


def eval_id(x):
    '''Define identity function for system of state n'''
    return x

def eval_gradient_id(x):

    return np.eye(len(x))

def eval_hessian_id(x):

    return np.zeros((len(x),len(x),len(x)))





        