# See http://matplotlib.org/1.4.3/examples/animation/moviewriter.html
# and https://stackoverflow.com/questions/16915966/using-matplotlib-animate-to-animate-a-contour-plot-in-python

from numpy import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import matplotlib.animation as manimation

sol_filename = "u.txt"
time_filename = "time.txt"
u = loadtxt(sol_filename)
t = loadtxt(time_filename)

print(f"max_(t,x,y) (u) = {amax(u)} (>1 => unstable)")

Nf = len(t)
Nx = u.shape[1]
Ny = int(u.shape[0]/Nf)

# print("N=%d" % N)
#u = u[:N, :]
#u.shape = (N,N)


x = linspace(-1, 1, Nx)
y = linspace(-1, 1, Ny)
j=1
X,Y = meshgrid(x, y)


fig = plt.figure()
ax = plt.axes(xlim=(-1, 1), ylim=(-1, 1), xlabel='x', ylabel='y')
# cvals = linspace(0,1,20)      # set contour values 
cont = plt.contourf(X, Y, u[:Ny,:])    # first image on screen
plt.colorbar()


def animate(i):
    global cont
    z = u[Ny*i : Ny*(i+1), :]
    for c in cont.collections:
        c.remove()  # removes only the contours, leaves the rest intact
    cont = plt.contourf(X, Y, z)
    plt.title(f't = {t[i]}')
    print(f"Frame {i}/{Nf} :: max(u) = {amax(z)}")
    return cont

anim = manimation.FuncAnimation(fig, animate, frames=Nf, repeat=False)
anim.save('animation.mp4', writer=manimation.FFMpegWriter(fps=15))