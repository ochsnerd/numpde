import mpl_toolkits.mplot3d
from numpy import *
from pylab import *

# u = loadtxt("u.txt")
# N = int(u.shape[1])
# print("N=%d" % N)
# u = u[:N, :]
# u.shape = (N,N)

sol_filename = "u.txt"
time_filename = "time.txt"
u = loadtxt(sol_filename)
t = loadtxt(time_filename)

Nf = len(t)
Nx = u.shape[1]
Ny = int(u.shape[0]/Nf)
x = linspace(-1, 1, Nx)
y = linspace(-1, 1, Ny)
j=1
X,Y = meshgrid(x, y)

fig = figure()
ax = fig.gca(projection="3d")
ax.plot_surface(X, Y, u[:Ny, :])
xlabel('x')
ylabel('y')
title('Initial condition')
show()