import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('./supercell_RDF.dat')

plt.figure(figsize = (8,6))

for i in range(100,1001,100):
    interval = data[data[:,0] == i][:,[1,2]]
    plt.plot(interval[:,0], interval[:,1],
             label = f'{i} ps', alpha = 0.5)

plt.title("Radial Distribution Function\n-From Random Tetrahedra-")
plt.xlabel("r ($10^{-10}$ m)")
plt.ylabel("RDF")

plt.legend()
plt.grid()

plt.show()
