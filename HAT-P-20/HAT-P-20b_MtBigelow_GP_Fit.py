""" HAT-P-20b: Mont4K 2018 Light Curve Fit (with CELERITE GP)
           ~ Arsh R. Nadkarni (UArizona), 2021 """  

#### Import Libraries ####
import numpy as np
from astropy.io import ascii
from astropy import units as u
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib import rcParams; rcParams["figure.dpi"] = 150
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import (AutoMinorLocator)
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
from astropy.table import Table, Column
from scipy.stats import norm
import emcee
import batman
import corner
import celerite
from celerite import terms

#### SETUP ####

# Import HAT-P-20b_Measurements.txt file and extract the light curve data.
data = ascii.read('/home/anadkarni/NRC_33-LW_Grism_Stability/HAT-P-20/HAT-P-20b_Measurements.txt')

# Extract Time, Flux and Error in Flux
time = data['JD_UTC']
flux = data['rel_flux_T1_dn']
flux_err = data['rel_flux_err_T1_dn']*3.124543156214643

#### BATMAN MODEL ####

# Initial Transit Parameters

t0 = 2455917.64431
t0_err = 0.00006
i = 86.3
i_err = 0.1
delta_F = 0.0173
delta_F_err = 0.0002
rp = np.sqrt(delta_F)
rp_err = 0.5*(delta_F_err/delta_F)
a_AU = 0.03671*u.au
a_AU_err = 0.00027*u.au
a_Rsun = a_AU.to(u.solRad)
a_Rsun_err = a_AU_err.to(u.solRad)
Rstar = 0.744*u.solRad
Rstar_err = 0.011*u.solRad
a_final = a_Rsun/Rstar
a_final_err = np.sqrt((a_Rsun_err/a_Rsun)**2 + (Rstar_err/Rstar)**2)
P = 2.8753172
P_err = 0.0000003

pars = np.zeros(7)
pars[0] = t0                       # t0
pars[1] = np.cos(i*(np.pi/180))    # cos(i)                        
pars[2] = rp                       # Rp/Rstar                    
pars[3] = np.log10(a_final)        # log(a/Rstar)
pars[4] = np.log10(P)              # log(P)
pars[5] = -4                       # log(a)[HP]
pars[6] = -13                      # log(c)[HP]

""" TRANSIT MODEL FUNCTION """
def transitmodel(pars, time):
    params = batman.TransitParams()
    params.t0 = pars[0]                            # t0
    params.inc = np.degrees(np.arccos(pars[1]))    # cos(i)
    params.rp = pars[2]                            # Rp/Rstar
    params.a = 10**(pars[3])                       # log(a/Rstar)
    params.per = 10**(pars[4])                     # log(P)
    params.ecc = 0.0136                            # e
    params.w = 90                                  # w                 
    params.limb_dark = "quadratic"       # Limb Darkening Model   
    params.u = [0.77822530, 0.019420697] # Limb Darkening Coefficients - EXOFAST
    transitmod = batman.TransitModel(params, time)
    lc = transitmod.light_curve(params)
    return lc

""" GP MODEL FUNCTION """
def gpmodel(pars,x,yerr,residuals):
    kernel = terms.RealTerm(log_a = pars[5], log_c = pars[6])
    gp = celerite.GP(kernel, mean=0.0)
    gp.compute(x, yerr)
    lnprobmodel = gp.log_likelihood(residuals)
    pred_mean, pred_var = gp.predict(residuals, x, return_var=True)
    return lnprobmodel, pred_mean, pred_var

#### MODEL CHECK ####

# Required Arrays 
model_flux_time = transitmodel(pars,time)
residual_flux_time = flux - model_flux_time
gpmodel_flux_time = gpmodel(pars,time,flux_err,residual_flux_time)
pm_time = gpmodel_flux_time[1]
time_init = (time-(pars[0]-((10**pars[4])/4)))/(10**pars[4]) 
phase = np.mod(time_init,1)
model_phase = np.linspace(-0.25,0.75,len(time))
model_times = (model_phase*(10**pars[4])) + pars[0]
model_flux_times = transitmodel(pars, model_times)
residual_flux_times = flux - model_flux_times
gpmodel_flux_times = gpmodel(pars,model_times,flux_err,residual_flux_times)
pm_times = gpmodel_flux_times[1]

# Plot Rough Transit Model
f, a = plt.subplots(2)
a[0].plot(time, flux, ls='None', marker='.', ms=0.9, zorder=0, color='black', label ="Observed Flux")
a[0].plot(time, model_flux_time+pm_time, 'r-', label ="BATMAN Model Flux")
a[0].set_xlabel(r'BJD$\mathregular{_{TDB}}$')
a[0].set_title(f'HAT-P-20b: Intensity v/s Time', fontweight='bold')
a[0].legend(loc='best', prop={'size': 6.8})
a[1].plot(phase-0.25, flux, ls='None', marker='.', ms=0.9, zorder=0, color='black')
a[1].plot(model_phase, model_flux_times, 'r-')
a[1].set_xlim(-0.05,0.05)
a[1].set_xlabel('Phase')
a[1].set_title(f'HAT-P-20b: Intensity v/s Phase', fontweight='bold')
for ax in a.flat:
    ax.set_ylabel('Intensity')
    ax.xaxis.set_minor_locator(AutoMinorLocator()) 
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='minor', length=2.5, color='k')
f.tight_layout(pad=1)
plt.savefig(f'premcmc/HAT-P-20b_PreMCMC+GP_LC.png')

#### MCMC ####

""" GELMAN-RUBIN STATISTIC (MCMC CONVERGENCE) FUNCTION """
def GelmanRubin(chains):
   nwalker = chains.shape[0]
   niter = chains.shape[1]
   npar = chains.shape[2]
   grarray = np.zeros(npar)
   for i in range(npar):
      sj2 = np.zeros(nwalker)
      chainmeans = np.zeros(nwalker)
      for j in range(nwalker):
         chainmeans[j] = np.mean(chains[j,:,i])
         sj2[j] = np.sum((chains[j,:,i]-chainmeans[j])**2.) / (niter-1)
      W = np.sum(sj2) / nwalker
      ThetaDoubleBar = np.sum(chainmeans) / nwalker
      B = np.sum((chainmeans-ThetaDoubleBar)**2.) * niter / (nwalker-1)
      VarTheta = (1-(1/niter))*W + (B/niter)
      grarray[i] = np.sqrt(VarTheta/W)
   return grarray

""" MCMC FUNCTION """
def RunMCMC(p0,nburn,nprod,nwalkers,ndim,log_probability,args):
    # Gather Results
    def medval_params(sampler,i): # args = sampler.flatchain[:,0/1/2...]
        samples = sampler.flatchain
        median = np.nanmedian(samples[:,i])
        minus = median - np.percentile(samples[:,i],16.0)
        plus = np.percentile(samples[:,i],84.0)- median
        return median, plus, minus
    from multiprocessing import Pool
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=args, pool=pool)
        print("Running burn-in...")
        p0,_,_ = sampler.run_mcmc(p0, nburn, progress=True)
        sampler.reset()
        print("Running production...")
        p0,_,_ = sampler.run_mcmc(p0, nprod, progress=True)
    np.savez(f'HAT-P-20b_MCMC+GP_Results', chain=sampler.chain, flatchain=sampler.flatchain, lnprob=sampler.lnprobability, flatlnprob=sampler.flatlnprobability)
    # Gelman-Rubin Statistic
    chains = sampler.chain
    gr = GelmanRubin(chains)
    print(f"Gelman-Rubin Statistic = {gr}")
    # Walker Performance
    f, a = plt.subplots(7, sharex=True)
    f.suptitle(f"HAT-P-20b: MCMC - Walker Performance", fontweight='bold')
    samples = sampler.get_chain()
    labels = ['$t_0$','cos(i)','$R_{p}/R_{*}$','$log(a/R_{*})$','log(P)','log(a)','log(c)'] #,'log(a)','log(c)','$\sqrt{e}sin(\omega)}$','$\sqrt{e}cos(\omega)}$','$M_{p}/M_{*}$','$\gamma$','$dv/dt$','$T_{eff}$','$R_{*}$','Parallax','A-V']
    for i in range(7):
        ax = a[i]
        ax.plot(samples[:, :, i], "k", alpha=0.3)
        ax.set_xlim(0, len(samples))
        ax.set_ylabel(labels[i], fontsize='x-small') 
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='minor', length=2.5, color='k')
    a[-1].set_xlabel("Step Number")
    plt.savefig(f"mcmc/HAT-P-20b_MCMC+GP_WP.png")
    # Fit Parameters
    pars_fit = np.zeros(7)
    pars_fit[0] = medval_params(sampler,0)[0]        # t0
    pars_fit[1] = medval_params(sampler,1)[0]        # cos(i)                  
    pars_fit[2] = medval_params(sampler,2)[0]        # Rp/Rstar                    
    pars_fit[3] = medval_params(sampler,3)[0]        # log(a/Rstar)
    pars_fit[4] = medval_params(sampler,4)[0]        # log(P)
    pars_fit[5] = medval_params(sampler,5)[0]        # log(a)
    pars_fit[6] = medval_params(sampler,6)[0]        # log(c)

    # Uncertainties on Parameters
    pars_fit_err = np.zeros(7)
    pars_fit_err[0] = (medval_params(sampler,0)[1]+ medval_params(sampler,0)[2])/2        # t0
    pars_fit_err[1] = (medval_params(sampler,1)[1]+ medval_params(sampler,1)[2])/2        # cos(i)                  
    pars_fit_err[2] = (medval_params(sampler,2)[1]+ medval_params(sampler,2)[2])/2        # Rp/Rstar                    
    pars_fit_err[3] = (medval_params(sampler,3)[1]+ medval_params(sampler,3)[2])/2        # log(a/Rstar)
    pars_fit_err[4] = (medval_params(sampler,4)[1]+ medval_params(sampler,4)[2])/2        # log(P)
    pars_fit_err[5] = (medval_params(sampler,5)[1]+ medval_params(sampler,5)[2])/2        # log(a)
    pars_fit_err[6] = (medval_params(sampler,6)[1]+ medval_params(sampler,6)[2])/2        # log(c)
    
    # Corner Plot
    flatchain = sampler.flatchain
    figure = corner.corner(flatchain, labels = ['$t_0$','cos(i)','$R_{p}/R_{*}$','$log(a/R_{*})$','log(P)','log(a)','log(c)'], label_kwargs=dict(fontsize=20), truths=[pars_fit[0], pars_fit[1], pars_fit[2], pars_fit[3], pars_fit[4], pars_fit[5], pars_fit[6]], truth_color='#ff0000', quantiles=[0.16, 0.5, 0.84], top_ticks=True, verbose=False)
    plt.savefig(f"mcmc/HAT-P-20b_MCMC+GP_Corner.png")
    # Return Required Values
    return pars_fit, pars_fit_err

# Define Priors
cosi_err = i_err*np.sin(i*(np.pi/180))
loga_err = a_final_err/(a_final*np.log(10))
logP_err = P_err/(P*np.log(10))

priors = np.zeros((2,7))
priors[:,0] = [t0, 2*t0_err]                         # t0
priors[:,1] = [np.cos(i*(np.pi/180)), cosi_err]      # cos(i)
priors[:,2] = [rp,rp_err]                            # Rp/Rstar
priors[:,3] = [np.log10(a_final),loga_err]           # log(a/Rstar)
priors[:,4] = [np.log10(P),2*logP_err]               # log(P)
priors_to_apply = [0,1,2,3,4]
print(priors)

######### Define log_probability Function #########
def log_probability(pars, priors, time, flux, flux_err, priors_to_apply):
    # Conditions
    if pars[1] < 0 or pars[1] > 1: return -np.inf
    if pars[5] > -3.0: return -np.inf
    if pars[6] > -12.0: return -np.inf
    # Calculations
    model = transitmodel(pars,time)
    residuals = flux - model
    kernel = terms.RealTerm(log_a = pars[5], log_c = pars[6])
    gp = celerite.GP(kernel, mean=0.0)
    gp.compute(time,flux_err)
    lnprobmodel_init = gp.log_likelihood(residuals)
    #lnprobmodel_init = -0.5*np.sum(((flux-model)**2/flux_err**2) + np.log(2*np.pi*flux_err**2))
    lnprobprior = 0
    for i in priors_to_apply:
        lnprobprior += -(pars[i] - priors[0][i])**2/(2*priors[1][i]**2) - np.log(np.sqrt(2*priors[1][i]**2*np.pi))
    lnprobmodel_total = lnprobmodel_init + lnprobprior
    return lnprobmodel_total

#### Run MCMC ####
args = (priors, time, flux, flux_err, priors_to_apply)
nwalkers = 14
ndim = 7
nburn = 20000 # No of Burn-In Steps
nprod = 40000 # No of Production Steps
scales = np.zeros(7)
scales[0] = 0.0008       # t0 
scales[1] = 0.008        # cos(i)                  
scales[2] = 0.001        # Rp/Rstar                    
scales[3] = 0.017        # log(a/Rstar)
scales[4] = 0.02         # log(P)
scales[5] = 0.1          # log(a)
scales[6] = 0.1          # log(c)
p0 = emcee.utils.sample_ball(pars, scales, nwalkers)
mcmc_init = RunMCMC(p0,nburn,nprod,nwalkers,ndim,log_probability,args)
pars_fit, pars_fit_err = mcmc_init