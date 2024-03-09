library(reticulate)
seed = 10 # the same in python
set.seed(seed)

# set path for NS-SG-sample.py
source_python('PDE_solver/NS-SG-sample.py')
source('Functions/all_functions.R')

num_xy = 5000; num_t = 60

snr_db_seq <- c(seq(20, 40, 2), Inf)
eta <- 10^(-snr_db_seq/20)
# simulate data
s_time <- Sys.time()
# Rudy gave 5000*30 we sample 5000*60
# system.time(RD_noise_data <- lapply(0:5*0.005, function(x){sam_rd(num_xy=5000,num_t=60,noise=x,cores=20)})) 
system.time(RD_noise_data <- lapply(eta/2, function(x){sam_rd(num_xy=5000,num_t=60,noise=x,cores=60)})) 
save(RD_noise_data,file=sprintf('Data/RD_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
e_time <- Sys.time()
difftime(e_time,s_time)
# Time difference of 19.02622 hours
names(RD_noise_data) <- sapply(snr_db_seq, function(x){paste0('SNR_dB=',x)})
save(RD_noise_data,file=sprintf('Data/RD_noise_data_%s*%s_seed_%s_snr.RData', num_xy, num_t, seed))
