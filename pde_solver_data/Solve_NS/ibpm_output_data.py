import numpy as np
import os

data_path = '/nmt/d/GitHub/ARGOS-Ral/pde_solver_data/Solve_NS/' # change this to your .plt data path

filenames = sorted([os.path.join(data_path,f) for f in os.listdir(data_path) if f[-3:] == 'plt'])
timesteps = len(filenames)

U = np.zeros((449,199,timesteps))
V = np.zeros((449,199,timesteps))
W = np.zeros((449,199,timesteps))

# save output files to .npy. Each of file has about 102 MB
np.save('U1.npy', U)
np.save('V1.npy', V)
np.save('W1.npy', W)