
import numpy as np
import matplotlib.pyplot as plt

T_gen = np.genfromtxt("heatmap.csv")
X = np.genfromtxt("coords.csv")
T1 = 300
T2 = 400

k1 = 45
k2 = 94

L = 1

T_gen = T_gen[0,:]
X = X[:T_gen.shape[0],0]

print(T_gen.shape,X.shape)


def T(x):
    if x < 0:
        return T1 - 2*(T1-T2)/(L*k1*(1/k1+1/k2))*(x+L/2)
    return T2 - 2*(T1-T2)/(L*k2*(1/k1+1/k2))*(x-L/2)

T_vec = np.vectorize(T)
Y = T_vec(X)


plt.plot(X,Y,c= 'r')
plt.plot(X,T_gen)
plt.xlabel("x [m]")
plt.ylabel("T [K]")
plt.legend(["analitical","FEM"])

plt.show()