# Creating PDE data 
This folder contains how to solve PDEs used in the paper.

## Most PDEs
Most simple PDEs are solved using **pde_data** function in *ode_int.py*. Some examples are shown in the *.py* file. We use this to solve Advection-diffusion, Burgers, Cable, Duffusion, Heat and Transport equations. 

## Burgers equation
We can use both *ode_int.py* and *burgers.m* to solve this equation. *burgers.mat* is created by *burgers.m*.

## KdV equation
We follow the book *Handbook of Nonlinear Partial Differential Equations* [1] (page 858) to solve this equation, codes are in *kdv_solver.R*. *kdv.mat* is created by Matlab, but uses the same codes in *kdv_solver.R*. 

## Schrodinger equation
The code is shown in *qho_solver.R* using a spectral method.

## Navier-Stokes & Reaction-diffusion
We follow the way that Rudy [2] used in PDE-FIND ([GitHub link](https://github.com/snagcliffs/PDE-FIND/tree/master/Datasets)) to solve both PDEs. 

*NS-SG-sample.py* is used to create the candidate library for Both Navier-Stokes and Reaction-diffusion. 

### Navier-Stokes 
The Navier-Stokes equation is solved using the Immersed Boundary Projection Method (IBPM) [3,4]. Please refer to the [GitHub repository](https://github.com/cwrowley/ibpm) for more information. All required files are located in the folder [*Navier-Stokes*](./Navier-Stokes/).
- [*NS_SG_data.R*](../Data/Navier-Stokes/NS_SG_data.R) runs Python functions in *NS-SG-sample.py* to create the candidate library and save it as *NS_noise_data_5500\*65_seed_10_snr.RData*, which occupies 2.13 GB. Running this *.Rdata* file, we used 22 CPU cores and each uses 2 threads with about 3.5 hours. The CPU is 2nd Gen AMD EPYC 7702.

### Reaction-diffusion
Reaction-diffusion is solved using *RD_solver.m*. The defined function is in *reaction_diffusion_rhs.m*. It may take about 6 hours to get the reaction-diffusion data, and the size output file *reaction_diffusion_big.mat* is about 767 MB. 
- [*RD_SG_data.R*](../Data/Reaction-diffusion/RD_SG_data.R) runs Python functions in *NS-SG-sample.py* to create the candidate library and save it as *RD_noise_data_5000\*60_seed_10_snr.RData*, which occupies 1.81 GB. Running this *.Rdata* file, we used 60 threads with about 31.5 hours. The CPU is 2nd Gen AMD EPYC 7702.

### References
[1] Andrei D. Polyanin and Valentin F. Zaitsev, Handbook of Nonlinear Partial Differential Equations, 2nd ed. New York: Chapman and Hall/CRC, 2012.

[2] S. H. Rudy, S. L. Brunton, J. L. Proctor, and J. N. Kutz, ‘Data-driven discovery of partial differential equations’, Sci. Adv., vol. 3, no. 4, p. e1602614, Apr. 2017, doi: 10.1126/sciadv.1602614.

[3] K. Taira and T. Colonius, ‘The immersed boundary method: A projection approach’, Journal of Computational Physics, vol. 225, no. 2, pp. 2118–2137, 2007, doi: 10.1016/j.jcp.2007.03.005.

[4] T. Colonius and K. Taira, ‘A fast immersed boundary method using a nullspace approach and multi-domain far-field boundary conditions’, Computer Methods in Applied Mechanics and Engineering, vol. 197, no. 25, pp. 2131–2146, 2008, doi: 10.1016/j.cma.2007.08.014.

