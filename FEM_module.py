import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
from shapes import shapes


class FemSolver():
    #geometry
    len_x          = 1
    len_y          = 1
    n_element_x          = 20
    n_element_y          = 10
    nodes_in_element      = 4
    
    def __init__(self):
        self.update()
        self.material = np.ones((self.n_element_y,self.n_element_x))

        self.dirchlet_conditions = np.zeros((self.n_element_y,self.n_element_x))

        self.neumann_conditions = np.zeros((self.n_element_y,self.n_element_x))

        self.dirchlet_conditions[0,:] = 10
        self.dirchlet_conditions[-1,:] = 100

    def set_materials(self,material_new):
        self.material = material_new
    
    def set_dirchlet(self, new_d):
        self.dirchlet_conditions = new_d 
    
    def set_neumann(self, new_d):
        self.neumann_conditions = new_d 


    def update(self):
        self.dx          = self.len_x/(self.n_element_x-1)
        self.dy          = self.len_y/(self.n_element_y-1)
        

    
        
        self.nex         = self.n_element_x-1
        self.ney         = self.n_element_y-1
        self.n_nodes        = self.n_element_x*self.n_element_y
        self.n_elements         = self.nex*self.ney
        self.global_coords      = np.zeros((self.n_nodes,2))

        self.T           = np.zeros(self.n_nodes) #initial T, not strictly needed


        
    def solve(self):
        
        
        self.update()

        if self.material.shape != (self.n_element_y,self.n_element_x):
            self.material = np.ones((self.n_element_y,self.n_element_x))

            self.dirchlet_conditions = np.zeros((self.n_element_y,self.n_element_x))

            self.neumann_conditions = np.zeros((self.n_element_y,self.n_element_x))

        #Setting global coordinates
        id = 0
        for i in range(0,self.n_element_y):
            for j in range(0,self.n_element_x):
                self.global_coords[id,0] = -self.len_x/2 + j*self.dx
                self.global_coords[id,1] = -self.len_y/2 + i*self.dy
                id          = id + 1

        # FEM connectivity
        # stores the node numbers belonging to each element
        elements_to_nodes   = np.zeros((self.n_elements,self.nodes_in_element), dtype=int)

        for i_element in range(0,self.n_elements):
            row        = i_element//self.nex   
            ind        = i_element + row
            elements_to_nodes[i_element,:] = [ind, ind+1, ind+self.n_element_x+1, ind+self.n_element_x]
            
        # Gauss integration points
        n_integration_points   = 4
        gauss = np.array([[ -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3), -np.sqrt(1/3)], [-np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3)]]).T.copy()

        # Storage
        Rhs_all = np.zeros(self.n_nodes)

        #indicies for the sparse matrix
        I       = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))
        J       = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))
        
        #global stiffnes matrix

        K       = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))

        for i_element in range(0,self.n_elements):

            element_coords = np.take(self.global_coords, elements_to_nodes[i_element,:], axis=0 )
            
            K_local    = np.zeros((self.nodes_in_element,self.nodes_in_element))

            #empty for the sake of completness
            Rhs_el = np.zeros(self.nodes_in_element)
            
            for i_integration_point in range(0,n_integration_points):        
                # 1. update shape functions
                xi      = gauss[i_integration_point,0]
                eta     = gauss[i_integration_point,1]
                #derivatives of the shape functions
                N, dNds = shapes(xi, eta)
                
                # 2. set up Jacobian, inverse of Jacobian, and determinant
                Jac     = np.matmul(dNds,element_coords) #[2,nnodel]*[nnodel,2]
                invJ    = np.linalg.inv(Jac)     
                detJ    = np.linalg.det(Jac)
                
                # 3. get global derivatives
                dNdx    = np.matmul(invJ, dNds) # [2,2]*[2,nnodel]
                
                #setting material properties
                grid_x = int ((np.mean(element_coords[:,0])/self.len_x+0.5) * self.n_element_x)
                grid_y = int ((np.mean(element_coords[:,1])/self.len_y+0.5) * self.n_element_y)

                k_element = self.material[grid_y,grid_x]

                K_local     = K_local + np.matmul(dNdx.T, dNdx)*detJ*k_element # [nnodel,1]*[1,nnodel] / weights are missing, they are 1
                
                # 5. assemble right-hand side, no source terms, just here for completeness
                Rhs_el     = Rhs_el + np.zeros(self.nodes_in_element)
            
            # assemble csr ids
            I[i_element,:]  =  (elements_to_nodes[i_element,:]*np.ones((self.nodes_in_element,1), dtype=int)).T.reshape(self.nodes_in_element*self.nodes_in_element)
            J[i_element,:]  =  (elements_to_nodes[i_element,:]*np.ones((self.nodes_in_element,1), dtype=int)).reshape(self.nodes_in_element*self.nodes_in_element)
            
            K[i_element,:]  =  K_local.reshape(self.nodes_in_element*self.nodes_in_element)
            
            Rhs_all[elements_to_nodes[i_element,:]] += Rhs_el

        A_all = csr_matrix((K.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element),(I.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element),J.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element))),shape=(self.n_nodes,self.n_nodes))


   

        #Dirchlet Conditions

        ind_dc_2d = np.where(self.dirchlet_conditions !=0)
        ind_dc = ind_dc_2d[1] + ind_dc_2d[0]*self.n_element_x
        val_dc = np.array([self.dirchlet_conditions[ind_dc_2d[0][i],ind_dc_2d[1][i]] for i in range(len(ind_dc_2d[0]))])

        Free    = np.arange(0,self.n_nodes)
        Free    = np.delete(Free, ind_dc)

        tmp_dc     = A_all[:,ind_dc]
        Rhs_all = Rhs_all - tmp_dc.dot(val_dc)

        Rhs_all +=  self.neumann_conditions.flatten()

        # solve reduced system
        self.T[Free] = spsolve(A_all[np.ix_(Free, Free)],Rhs_all[Free])
        self.T[ind_dc] = val_dc

        return self.T.reshape((self.n_element_y,self.n_element_x))


if __name__ == "__main__":
    fs = FemSolver()

    T = fs.solve()

    plt.matshow(T)
    plt.show()