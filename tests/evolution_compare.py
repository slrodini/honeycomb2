import numpy as np

root_folder = "/mnt/d/Dropbox/ThinkPad/Documents/Code/honeycomb2/"

data_old = np.loadtxt(root_folder + "dat/T_and_DT_final_scale_1p0e+4.dat")

data_new = np.loadtxt(root_folder + "build/T_and_DT_final_scale.dat")

m = 0

for i in range(data_new.shape[0]):
    for j in range(data_new.shape[1]):
        m = np.amax([m, np.abs(data_new[i, j] - data_old[i, j])])
print("Max difference: ", m)
