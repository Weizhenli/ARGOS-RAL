import numpy as np
import os
import time

data_path = '~/GitHub/ARGOS-RAL/Data/Navier-Stokes/plt_files/'

filenames = sorted([os.path.join(data_path,f) for f in os.listdir(data_path) if f[-3:] == 'plt']) # find all .plt files
timesteps = len(filenames)

U = np.zeros((449,199,timesteps))
V = np.zeros((449,199,timesteps))
W = np.zeros((449,199,timesteps))

for timestep in range(timesteps):
    timestep_data = np.genfromtxt(filenames[timestep], delimiter=' ',skip_header=6)
    for i in range(449):
        for j in range(199):
            U[i,j,timestep] = timestep_data[i+449*j, 2]
            V[i,j,timestep] = timestep_data[i+449*j, 3]
            W[i,j,timestep] = timestep_data[i+449*j, 4]

# save as npy
np.save('U1.npy', U)
np.save('V1.npy', V)
np.save('W1.npy', W)
