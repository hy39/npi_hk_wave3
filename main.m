% Main function to plot the figures:
filename = 'out/mcmc/m11.1/mcmc_output_m11.1(1).mat';
%filename = 'out/mcmc/m11.1/mcmc_output_m11.1(2).mat'; %MCMC second
%tracjectory
dat = open('filename);
plot_HK_dynamics_soc11_1_sliding(dat.PosteriorSamples, dat.sys_par, dat.par);