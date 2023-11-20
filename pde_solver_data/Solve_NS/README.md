# Solve Navier-Stokes equation
Please download or git the [IBPM](https://github.com/cwrowley/ibpm) to this folder. 

- *ibmp.geom* defines the geometry of the cylinder.
- *ibmp.sh* is a shell file for submitting a task in Linux.
- *ibmp15300.txt* is the output file from IBPM at the timestep 15300. The original output file ends with *.plt*, but we changed it to *.txt* to allow reading into Python and R.
- *ibmp.txt* is the command line for running IBPM. 
- *ibpm_output_data.py* transform *.plt* data to *.npy*, output *U1.npy*, *V1.npy* and *W1.npy*. These *.npy* files are used in [*NS_SG-sample.py*](../NS-SG-sample.py) to create the candidate library. 