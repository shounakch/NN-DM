#setwd("~/NNDP/NNDP Multivariate Codes/Updated Codes/New_Simulations_NonDisjoint_NNDP/n=200_p=2_dnum=1")

## Summary NNDP ##
n = 200
dnum = 20
p = 6
dir.name = paste("n=", n, "_p=", p, "_dnum=", dnum, sep = "")
wd = paste("~/NNDP/NNDP Multivariate Codes/Updated Codes/New_Simulations_NonDisjoint_NNDP/", dir.name, sep = "")
setwd(wd)
wd
av = (1:20) ##Has to be supplied
n_replicates = length(av)
MC = 1000
l1_error = matrix(0, n_replicates, 2)
colnames(l1_error) = c("NNDP_median", "NNDP_mean")

i=1

for(i in 1:length(av))
{
  fname = paste("run", av[i], ".csv", sep = "")
  sim_data = read.csv(fname)
  sim_data = sim_data[,-1]
  f0 = sim_data[,(p+1)]
  #f0 = dberdev(sim_data[,1],dnum)
  f_NNDP_mean = rowMeans(sim_data[,(p+2):(p+MC+1)])
  f_NNDP_median = apply(sim_data[,(p+2):(p+MC+1)], 1, median)
  #f_DP = sim_data[,(p+MC+2)]
  #f_KDE = sim_data[,(p+MC+5)]
  df = data.frame(f0, f_NNDP_median, f_NNDP_mean)
  
  #l1_error = rep(0, 4)
  for(l in 1:2)
  {
    l1_error[i, l] = mean(abs(df[,1] - df[,(l+1)])/df[,1])
  }
}

l1_error

colMeans(l1_error)
