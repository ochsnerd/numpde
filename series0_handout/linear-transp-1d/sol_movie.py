# See http://matplotlib.org/1.4.3/examples/animation/moviewriter.html


from numpy import *
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import matplotlib.animation as manimation

try:
	sol_filename = sys.argv[1]
except:
	print "Usage: python sol_movie.py <name of the file containing data for u>"
	sys.exit()

u=loadtxt(sol_filename)
t = loadtxt("time" + sol_filename[1:])
x = linspace(0,5, u.shape[1])
Nf = u.shape[0]
j=1







FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title=sol_filename, artist='Burgers',
        comment='Look at the solution go!')
writer = FFMpegWriter(fps=15, metadata=metadata)

fig = plt.figure()
l, = plt.plot([], [], 'k-o')


plt.xlim([0,5])
plt.ylim([amin(u),amax(u)])
plt.xlabel('$x$')
plt.ylabel('$u(x)$')




with writer.saving(fig, sol_filename.replace('.txt','.mp4'), 100):
    for i in range(0, Nf, int(Nf/100)):
        plt.title("T = " + str(t[i]))
        l.set_data(x, u[i,:])
        writer.grab_frame()


