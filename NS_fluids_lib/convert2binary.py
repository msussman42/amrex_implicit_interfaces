import numpy as np

exact_f = np.loadtxt('exact_f.dat',dtype=np.float64)
exact_angle = np.loadtxt('exact_angle.dat',dtype=np.float64)
exact_centroid = np.loadtxt('exact_centroid.dat',dtype=np.float64)
initial_angle = np.loadtxt('initial_angle.dat',dtype=np.float64)

np.save('exact_f.npy',exact_f)
np.save('exact_angle.npy',exact_angle)
np.save('exact_centroid.npy',exact_centroid)
np.save('initial_angle.npy',initial_angle)