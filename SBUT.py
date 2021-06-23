# In[]
import matplotlib.pyplot as plt
import numpy as np
#from numpy.linalg import inv
import matplotlib
from sympy import Matrix, pprint
#matplotlib.use('PDF')
plt.rc('text', usetex=True)

# <md>
# Directed network with a master cluster and a slave cluster

# In[]
# Adjacency matrix A and diagonal matrices C1 and C2 encoding clusters
A = np.array([[0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 1, 0], [1, 0, 0, 1, 0, 0, 1], [0, 1, 0, 1, 0, 0, 1], [0, 0, 1, 0, 1, 1, 0]])
C1 = np.diag([1, 1, 1, 0, 0, 0, 0])
C2 = np.diag([0, 0, 0, 1, 1, 1, 1])

# In[]
# Convert adjacency matrix A into Laplacian matrix L
L = np.diag(np.sum(A, axis=1)) - A

# In[]
# Generate a single matrix M from L, C1 and C2
M = Matrix(2*C1@L@C1 + 5*C2@L@C2)

# In[]
# Find the basis for the common block upper triangular form by finding the Jordan decomposition of M
P, J = M.jordan_form()
TempL = P.inv()*L*P
pprint(J)

# In[]
# Reorder the basis P so that it reveals the common block upper triangular form
order = [0,3,4,5,6,1,2]
FinalL = TempL[order,order]
FinalP = P[:,order]#.simplify()
pprint(FinalL)
pprint(FinalP)


# <md>
# Directed network formed by four clusters with directional flow

# In[]
# Adjacency matrix A and diagonal matrices C1 to C4 encoding clusters
A = np.array([[0, 1, 0, 0, 0, 0, 0, 1], [1, 0, 1, 0, 0, 0, 1, 0], [1, 0, 0, 0, 0, 0, 0, 1], [1, 1, 0, 0, 1, 0, 0, 0], [0, 1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 1],[0, 0, 0, 0, 0, 1, 1, 0]])
C1 = np.diag([1, 1, 1, 0, 0, 0, 0, 0])
C2 = np.diag([0, 0, 0, 1, 1, 0, 0, 0])
C3 = np.diag([0, 0, 0, 0, 0, 1, 0, 0])
C4 = np.diag([0, 0, 0, 0, 0, 0, 1, 1])

# In[]
# Convert adjacency matrix A into Laplacian matrix L
L = np.diag(np.sum(A, axis=1)) - A

# In[]
# Generate a random combination M of all generating matrices so we don't have to worry about degeneracy (or use the form suggested in the paper for better alignment with clusters)
M = Matrix(2*C1@L@C1 + 3*C2@L@C2 + 5*C3@L@C3 + 7*C4@L@C4)

# In[]
# Find the basis for the common block upper triangular form by finding the Jordan decomposition of M
P, J = M.jordan_form()
TempL = P.inv()*L*P

# In[]
# Use the BlockSort algorithm to reorder the basis P so that it can be used to reveal the common block upper triangular form
order = [0,4,1,5,6,2,3,7]
FinalL = TempL[order,order]
FinalB = P[:,order].simplify()
pprint(FinalL)
print(FinalB)

# <md>
# Two intertwined clusters with bi-directional flow

# In[]
# Adjacency matrix A and diagonal matrices C1 and C2 encoding clusters
A = np.array([[0, 1, 0, 1, 0, 1], [1, 0, 1, 0, 0, 1], [0, 1, 0, 1, 1, 0], [1, 0, 1, 0, 1, 0], [1, 0, 0, 1, 0, 0], [0, 1, 1, 0, 0, 0]])
C1 = np.diag([1, 1, 1, 1, 0, 0])
C2 = np.diag([0, 0, 0, 0, 1, 1])

# In[]
# Convert adjacency matrix A into Laplacian matrix L
L = np.diag(np.sum(A, axis=1)) - A

# In[]
# Generate a random combination M of all generating matrices so we don't have to worry about degeneracy (or use the form suggested in the paper for better alignment with clusters)
#M = Matrix(2*C1@L@C1 + 5*C2@L@C2)
M = Matrix(2*L + 3*C1 + 5*C2)

# In[]
# Find the basis for the common block upper triangular form by finding the Jordan decomposition of M
P, J = M.jordan_form()
TempL = P.inv()*L*P

# In[]
# Use the BlockSort algorithm to reorder the basis P so that it can be used to reveal the common block upper triangular form
order = [4,5,3,0,1,2]
FinalL = TempL[order,order]
FinalB = P[:,order].simplify()
pprint(FinalL)
print(FinalB)
