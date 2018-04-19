import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
K = 1

def space_grid(a, b, NS):
    X = np.linspace(a, b, NS + 1)
    h = (b - a) / NS
    return X, h
def time_grid(a, b, NT):
    T = np.linspace(a, b, NT + 1)
    tau = (b-a) / NT
    return T, tau
def u_initial(x):
        return np.exp(-(x - 0.25) ** 2 / 0.01) + 0.1 * np.sin(20 * np.pi * x)
def f(x, t):
    return np.zeros(len(x))

X, h = space_grid(0, 1, 100)
T, tau = time_grid(0, 0.1, 10000)
r = K * tau / (h ** 2)
N = len(X)
M = len(T)
U = np.zeros((N, M))
U[:, 0] = u_initial(X)
U[0, :] = np.zeros(len(T))
U[-1, :] = np.zeros(len(T))

def forward(N, M):
    d = 1 - 2 * np.ones(N - 2) * r
    c = np.ones(N - 3) * r
    A = np.diag(c, -1) + np.diag(d) + np.diag(c, 1)
    for i in range(1, M):
        rhs = tau * f(X, T[i])
        rhs[1] = rhs[1] + r * U[0, i-1]
        rhs[-2] = rhs[-2] + r * U[-1, i-1]
        U[1: -1, i] = np.dot(A, U[1: -1, i - 1]) + rhs[1:-1]
    return U
U = forward(N, M)
U = U.T
X, T = np.meshgrid(X, T)
fig = plt.figure()
ax = Axes3D(fig)
ax.plot_surface(X, T, U, rstride=10, cstride=10, cmap = cm.viridis)
plt.xlabel('X')
plt.ylabel('T')
ax.set_zlabel('U(X, T)')
plt.show()