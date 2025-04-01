import numpy as np
import emcee
# import scipy.special as sc
# import scipy.stats as stats

class Mixture_Model():

    def __init__(self, **kwargs):
        self.ndim_check = kwargs.get('ndim_check',0)
        self.nwalkers = kwargs.get('nwalkers',20) # number of MCMC walkers
        self.nburn    = kwargs.get('nburn',1000)    # "burn-in" period to let chains stabilize
        self.nsteps   = kwargs.get('nsteps',4000) # number of MCMC steps to take
        self.nthreads = kwargs.get('nthreads',1)  # number of threads
        self.vlos_min = kwargs.get('vlos_min', -500)
        self.vlos_max = kwargs.get('vlos_max', 500)
        self.vlos_start = kwargs.get('vlos_start', 0)
        
        self.rs_min = kwargs.get('rs_min', -2)
        self.rs_max = kwargs.get('rs_max', 1)
        self.rhos_min = kwargs.get('rhos_min', 4)
        self.rhos_max = kwargs.get('rhos_max', 12)
        self.progress = kwargs.get('progress', True)
       
        self.rhalf = kwargs.get('rhalf', 1.)
        self.rhalf_error = kwargs.get('rhalf_error', 1.)
        self.ellipticity = kwargs.get('ellipticity', 0.)
        self.ellipticity_error = kwargs.get('ellipticity_error', .1)
        self.distance = kwargs.get('distance', 1.)
        self.distance_error = kwargs.get('distance_error', 1.)
        self.pmra = kwargs.get('pmra', 0.)
        self.pmra_error = kwargs.get('pmra_error', 0.)
        self.pmdec = kwargs.get('pmdec', 0.)
        self.pmdec_error = kwargs.get('pmdec_error', 0.)
        self.ra = kwargs.get('ra', 0.)
        self.dec = kwargs.get('dec', 0.)

        self.profile_option = kwargs.get('profile_option', 0)

        self.limit_distance_min = np.max((0., self.distance-7.*self.distance_error))
        self.limit_distance_max = self.distance+7.*self.distance_error

        self.limit_rhalf_min = np.max((0., self.rhalf-7.*self.rhalf_error))
        self.limit_rhalf_max = self.rhalf+7.*self.rhalf_error

        self.limit_ellipticity_min = 0.
        self.limit_ellipticity_max = 1.

        self.limit_pmra_min = self.pmra-7.*self.pmra_error
        self.limit_pmra_max = self.pmra+7.*self.pmra_error
        
        self.limit_pmdec_min = self.pmdec-7.*self.pmdec_error
        self.limit_pmdec_max = self.pmdec+7.*self.pmdec_error

        # self.create_starting_prior(self)
        self.data_vlos = kwargs.get('data_vlos')
        self.data_vlos_error = kwargs.get('data_vlos_error')
        self.data_ra = kwargs.get('data_ra')
        self.data_dec = kwargs.get('data_dec')


        self.starting_guesses, self.limits, self.parameters_name, self.ndim = self.create_starting_prior()
        print("self.ndim", self.ndim)
        self.sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim, self.lnprob_sigma_i,threads=self.nthreads)
        

    def create_starting_prior(self, ):
        ndim = 1
        sigma_ball, sigma_ball_vlos = 0.001, 1

        X = stats.truncnorm((self.vlos_min - self.vlos_start) / sigma_ball_vlos, (self.vlos_max - self.vlos_start) / sigma_ball_vlos, loc=self.vlos_start, scale=sigma_ball_vlos)
        start_vlos = X.rvs(self.nwalkers)
        starting_guesses = [ start_vlos]

        limits = [self.vlos_min, self.vlos_max]
        parameters_name=['vlos'] 

        if ndim != self.ndim_check:
            print("wrong ndim", ndim, self.ndim_check)
            
        return starting_guesses, limits, parameters_name, ndim
    
    def lnprior_sigma_i(self, theta,  ):
        return_prior = 0
        place=0

        for i in range(int(len(self.limits)/2)):
            if np.logical_or(theta[i] < self.limits[2*i], theta[i] > self.limits[2*i+1]):
                return -np.inf

        return return_prior

    def lnprob_sigma_i(self, theta  ):
        lp = self.lnprior_sigma_i(theta, )
        if not np.isfinite(lp):
            return -np.inf
        return lp + self.lnlike_sigma_i(theta  )
    
    def lnlike_sigma_i(self, theta, ):

        place=0
        mcmc_vlos = theta[place]
        place +=1

        if np.isnan(term ) or np.isinf(term):
            return -1e20

        return term
    
    def run_sampler(self, starting_guesses, nburn, nsteps):

        starting_guesses2 = np.vstack(self.starting_guesses).T

        print(starting_guesses2.shape)

        pos,prob,state = self.sampler.run_mcmc(starting_guesses2, nburn+ nsteps,progress=self.progress,) 

        return self.sampler