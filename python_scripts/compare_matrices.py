# in the following file we read two sparse matrices from a file and read them
import matplotlib.pyplot as pyplot
from numpy import*

x1,y1,z1 = loadtxt("../printed_matrices/global_matrix1",unpack=True);
x2,y2,z2 = loadtxt("../printed_matrices/global_matrix2",unpack=True);


diff1 = z1 - z2;


for i in range(0,len(diff1)):
	if(abs(diff1[i]) > 1e-5):
		print i;

