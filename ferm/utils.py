import numpy as np
from scipy import stats
from arspy.ars import adaptive_rejection_sampling


def sample_max_distribution(mu_x,sigma, pop, alpha):
    """ Compute boundaries and sample max distributions """
    mm = (np.sqrt(np.log(pop**2 / (2 * np.pi * np.log(pop**2 / (2 * np.pi))))) * (1 + 0.5772156649 / np.log(pop))) * sigma + mu_x
    ss = np.sqrt(sigma**2 * np.pi**2 / (12 * np.log(pop)))
    
    low = mm - alpha * ss
    up = mm + alpha * ss
    
    return low, up

def gaussian_distribution_max(sigma,mu,n):
    #print("n= ", n)
    #print("mu",mu)
    "Sample the maximum of n gaussian distributions with Adaptative rejection sampling"
    if n <=5: #in case of small population make it naive way to avoid problem of boundaries in the rejection algorithm
        return np.max(stats.norm.rvs(loc=mu,scale=sigma,size=n))
    else:
        alpha=3
        a,b = sample_max_distribution(mu,sigma, n, alpha)
        #print(a," ",b)
        domain = (float("-inf"), float("inf"))
        log_pdf_max_gaussian = lambda x,mu=mu,sigma=sigma,n=n: np.log((n/sigma)*stats.norm.pdf(x,loc=mu,scale=sigma)*(stats.norm.cdf(x,loc=mu,scale=sigma))**(n-1))
        ars_sample = adaptive_rejection_sampling(log_pdf_max_gaussian, a=a, b=b, domain=domain, n_samples=1)[0]
        return ars_sample
