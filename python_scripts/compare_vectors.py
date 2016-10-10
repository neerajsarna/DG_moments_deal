# in the following file we read two sparse matrices from a file and read them
import matplotlib.pyplot as pyplot
from numpy import*

x1 = loadtxt("../printed_matrices/system_rhs1",unpack=True);
x2 = loadtxt("../printed_matrices/system_rhs2",unpack=True);

diff1 = x1 - x2;

for i in range(0,len(diff1)):
	if(abs(diff1[i]) > 1e-5):
		print diff1[i];

