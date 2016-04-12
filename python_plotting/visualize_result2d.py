import matplotlib.pyplot as plt
from numpy import*
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import os

for fn in os.listdir('.'):
	if os.path.isfile(fn):
		if (("solution" in fn) & (not".eps" in fn)):
			print fn;
			data = loadtxt(fn,unpack=True);
			pnts_interpolate = 1000;

			# create points on the x and the y axis
			x = linspace(-1,1,pnts_interpolate);
			y = linspace(-1,1,pnts_interpolate);

			# points where the solution has to be interpolated
			grid_x,grid_y = meshgrid(x,y);

			# value of the solution at the interpolation points
			grid_solution = griddata(data[[0,1],:].transpose(),data[2,:].transpose(),(grid_x,grid_y),method='linear');

			CS=plt.contourf(grid_x,grid_y,grid_solution,100);
			cbar = plt.colorbar(CS);
			plt.xlabel('x');
			plt.ylabel('y');
			plt.title('u(x,y) for poisson equation');
			plt.show();

