from numpy import *
from pylab import *

import matplotlib.tri
from math import atan2

import sys

basename = sys.argv[1]
vertices = loadtxt("%s_vertices.txt" % basename)
indices = loadtxt("%s_triangles.txt" % basename)
values = loadtxt("%s_values.txt" % basename)

print("using %d vertices" % len(vertices))

grid = matplotlib.tri.Triangulation( vertices[:,0], vertices[:,1], indices)

tripcolor(grid, values);

colorbar()
show()
