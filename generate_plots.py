import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Plot MDS eigenvalues for some symmetric groups
for i in range(4,8):
    eigs = np.genfromtxt('S' + str(i) + '_Kendall_eigenvalues.csv', delimiter=',')
    plt.plot(eigs)
    plt.title("Eigenvalues for the MDS operator on " + r'$S_'+str(i)+'$') 
    plt.show()
