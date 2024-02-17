import numpy as np
from scipy.sparse.linalg import spsolve
from scipy.sparse import csr_matrix



class FemSolver():
    #geometry properties

    # length in x, y direction in m
    len_x = 1 
    len_y = 1

    # number of nodes in x, y
    n_node_x = 20
    n_node_y = 20

    # nodes per element, in this case 4 since we are using a rectangular
    nodes_in_element = 4
    
    def __init__(self):

        # arrays to store information about material, openings and initial conditions
        self.material = np.ones((self.n_node_y,self.n_node_x))

        self.dirchlet_conditions = np.zeros((self.n_node_y,self.n_node_x))

        self.neumann_conditions = np.zeros((self.n_node_y,self.n_node_x))


        # this specifies which nodes should be included in the calculation
        # every node that corresponds to a 0 is treated as non existing
        self.node_space = np.ones((self.n_node_y,self.n_node_x))


        # recalculating parameters that depend on others
        self.update()
        

        # this is just for demonstration, it gets overwritten anyway in the main app
        self.dirchlet_conditions[0,:] = 10
        self.dirchlet_conditions[-1,:] = 100


    def update(self):

        # length of individual elements
        self.dx = self.len_x/(self.n_node_x-1)
        self.dy = self.len_y/(self.n_node_y-1)
        

        self.n_elements_x = self.n_node_x-1
        self.n_elements_y = self.n_node_y-1


        # calculate the number of nodes
        self.n_nodes = np.count_nonzero(self.node_space)
        
        # create a matrix which specifies which elements are present and 
        # calculate the total number of elements
        self.element_space = self.node_space[:-1,:-1]*self.node_space[1:,:-1]*self.node_space[:-1,1:]*self.node_space[1:,1:]
        self.n_elements = np.count_nonzero(self.element_space)
        
        # initialize an array which will store the coordinates of all nodes for later
        self.global_coords = np.zeros((self.n_nodes,2))

        self.T = np.zeros(self.n_nodes) #initial T


        
    def solve(self):
        
        
        # check if the shapes of arrays match with the specified shape 
        # if not generate new empty arrays 
        if self.material.shape != (self.n_node_y,self.n_node_x):
            self.material = np.ones((self.n_node_y,self.n_node_x))

            self.dirchlet_conditions = np.zeros((self.n_node_y,self.n_node_x))

            self.neumann_conditions = np.zeros((self.n_node_y,self.n_node_x))
            self.node_space = np.ones((self.n_node_y,self.n_node_x))

        # recalculate all dependant variables
        self.update()
        

        #Setting up global coordinates

        # This array allows us to retrieve a node id based on its position

        coords_to_nodes = np.zeros((self.n_node_y,self.n_node_x))
        

        # this loop fills the array of global coordinates with coordinates
        # of subsequent nodes, if the node is marked as exisiting by the node_space array 
        id = 0
        for i in range(0,self.n_node_y):
            for j in range(0,self.n_node_x):
                
                if self.node_space[i,j] == 0:
                    continue
                self.global_coords[id,0] = -self.len_x/2 + j*self.dx
                self.global_coords[id,1] = -self.len_y/2 + i*self.dy

                coords_to_nodes[i,j] = id
                id          = id + 1

        # this will allow us to get the node numbers corresponding to each element
        elements_to_nodes   = np.zeros((self.n_elements,self.nodes_in_element), dtype=int)


        # iterate over elements and if the element exists add its nodes to the array
        id_element = 0
        for i in range(0,self.n_elements_y):
            for j in range(0,self.n_elements_x):

                if self.element_space[i,j] == 0:
                    continue

                node_1 = coords_to_nodes[i,j]
                node_2 = coords_to_nodes[i,j+1]
                node_3 = coords_to_nodes[i+1,j+1]
                node_4 = coords_to_nodes[i+1,j]

                elements_to_nodes[id_element,:] = [node_1,node_2,node_3,node_4]
                id_element += 1

        
        # points for the shape functions
        n_integration_points = 4
        gauss = np.array([[ -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3), -np.sqrt(1/3)], [-np.sqrt(1/3), -np.sqrt(1/3), np.sqrt(1/3), np.sqrt(1/3)]]).T.copy()

        # Generate empty matricies for later

        #indicies for the sparse matrix
        I = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))
        J = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))
        
        # matricies for the equation K * T = F
        K = np.zeros((self.n_elements,self.nodes_in_element*self.nodes_in_element))
        F = np.zeros(self.n_nodes)


        # iterate over all elements and for each perform all the necessary calculations 
        for i_element in range(0,self.n_elements):


            # retrieve coordinates of the nodes of the element
            element_coords = np.take(self.global_coords, elements_to_nodes[i_element,:], axis=0 )
            
            # generate the local stiffness matrix for later 
            K_local    = np.zeros((self.nodes_in_element,self.nodes_in_element))

  
            # iterate over each point and perform the necessary calculations
            for i_integration_point in range(0,n_integration_points):       


                # get values of local corrdinates xi and eta and the values of their derivatives
                xi      = gauss[i_integration_point,0]
                eta     = gauss[i_integration_point,1]
                N, dNds = self.shapes(xi, eta)
                
                # calculate the Jacobian its inverse and its determinant
                Jac     = np.matmul(dNds,element_coords) 
                invJ    = np.linalg.inv(Jac)     
                detJ    = np.linalg.det(Jac)
                
                # get global derivatives
                dNdx    = np.matmul(invJ, dNds)
                
                # calculate the node position and within the grid and get the corresponding value
                # of k from the material matrix
                grid_x = int ((np.mean(element_coords[:,0])/self.len_x+0.5) * self.n_node_x)
                grid_y = int ((np.mean(element_coords[:,1])/self.len_y+0.5) * self.n_node_y)

                k_element = self.material[grid_y,grid_x]

                # calculate the local stiffness matrix
                K_local     = K_local + np.matmul(dNdx.T, dNdx)*detJ*k_element 
                
         
            
            # assemble csr ids
            I[i_element,:]  =  (elements_to_nodes[i_element,:]*np.ones((self.nodes_in_element,1), dtype=int)).T.reshape(self.nodes_in_element*self.nodes_in_element)
            J[i_element,:]  =  (elements_to_nodes[i_element,:]*np.ones((self.nodes_in_element,1), dtype=int)).reshape(self.nodes_in_element*self.nodes_in_element)
            
            K[i_element,:]  =  K_local.reshape(self.nodes_in_element*self.nodes_in_element)
            
        # convert the global stiffness matrix into a sparse matrix
        K_csr = csr_matrix((K.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element),(I.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element),J.reshape(self.n_elements*self.nodes_in_element*self.nodes_in_element))),shape=(self.n_nodes,self.n_nodes))


   

        # Dirchlet Conditions

        # Find which elements are marked to be set by the dirchlet conditions
        # and generate an array of corresponding temperatures
        ind_dc_2d = np.where(self.dirchlet_conditions !=0)
        ind_dc = coords_to_nodes[ind_dc_2d].astype(int)
        val_dc = np.array([self.dirchlet_conditions[ind_dc_2d[0][i],ind_dc_2d[1][i]] for i in range(len(ind_dc_2d[0]))])


        # indicies marked by this will be excluded from the final calculation
        # since we already know their values 
        bounded = np.arange(0,self.n_nodes)
        bounded = np.delete(bounded, ind_dc)

        # adjust the F matrix for the set temperatures
        tmp_dc = K_csr[:,ind_dc]
        F = F - tmp_dc.dot(val_dc)

        # add the Neumann conditions the the F matrix
        non_zero = np.where(self.node_space.flatten()!= 0)
        F +=  self.neumann_conditions.flatten()[non_zero]

        # solve the system while ingnoring the nodes which are set by the Dirchlet conditions
        self.T[bounded] = spsolve(K_csr[np.ix_(bounded, bounded)],F[bounded])
        self.T[ind_dc] = val_dc

        # reshape the flattened array back into a matrix and fill the nodes which are marked as 
        # non exisitant with zeroes
        T_reshaped = np.zeros((self.n_node_y,self.n_node_x))
        id = 0
        for i in range(0,self.n_node_y):
            for j in range(0,self.n_node_x):
                if self.node_space[i,j] == 0:
                    continue
                T_reshaped[i,j] = self.T[id]
                id += 1

        return T_reshaped
    
    def shapes(self,xi, eta):
    
        #shape functions and their derivatives
        
        N1 = 0.25*(1-xi)*(1-eta)
        N2 = 0.25*(1+xi)*(1-eta)
        N3 = 0.25*(1+xi)*(1+eta)
        N4 = 0.25*(1-xi)*(1+eta)

        N = np.array([N1, N2, N3, N4])
        
        dNds = np.zeros((2,4))
    
        dNds[0,0]   =  0.25*(-1  + eta) #derivative with xi
        dNds[1,0]   =  0.25*(xi  -  1) #derivative with eta

        #derivatives of second shape function with local coordinates
        dNds[0,1]   =  0.25*(1   - eta)
        dNds[1,1]   =  0.25*(-xi -  1)

        #derivatives of third shape function with local coordinates
        dNds[0,2]   =  0.25*(eta  +  1)
        dNds[1,2]   =  0.25*(xi  +  1)

        #derivatives of fourth shape function with local coordinates
        dNds[0,3]   =  0.25*(-eta -  1)
        dNds[1,3]   =  0.25*(1   - xi)
        
        return N, dNds

if __name__ == "__main__":
    # this can be used in case the solver is to be used independently
    pass