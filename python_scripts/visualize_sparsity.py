import matplotlib.pyplot as plt 
from numpy import*

x,y = loadtxt("sparsity_pattern1",unpack=True);

plt.plot(x,y,"o");
plt.title("sparsity pattern");
plt.grid(True);
plt.show();
