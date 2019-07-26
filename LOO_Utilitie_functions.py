#########################           IMPORTS         ###################################
import numpy as np
import matplotlib.pyplot as plt
import pystan as ps
import sys
sys.path.append("STAN_Utilities/")
import stan_utility
import psis
from scipy.special import logsumexp
import statistics as st
import scipy

#########################           UTILITY FUNCTIONS        ###################################

def eqrate(raterange,t) :
    """ Convert a list of earthquakes time of occurence (t) to a list of centered interval time (raterange is its value in days) + a list of rates over the interval time. """

    length=raterange
    nbint= int(t[-1]//length)
    inttime=[length/2+i*length for i in range(0,nbint)] # list of times
    intrate=[] # list of rates
    for j in range(0,nbint) : # itering the list t
        start=j*length # start time value
        end=(j+1)*length # end time value
        ind=[]
        for i in range(len(t)) :
            if start<t[i]<end :
                ind.append(i) # earthquakes occuring within the time interval
        intrate.append(len(ind))

    if t[-1]%length!=0 : # if the total time is not divisble by the time interval
        # add a last smaller interval
        inttime.append(nbint*length+(t[-1]-nbint*length)/2)
        start=nbint*length # start time value
        end=t[-1] # end time value
        ind=[]
        for i in range(len(t)) :
            if start<t[i]<end :
                ind.append(i) # earthquakes occuring within the time interval
        intrate.append(len(ind))

    return inttime,intrate



def cpsearch (count_data,R) :
    """ Return the probability density function of the change point over the list of rates count_data"""

    N=len(count_data)#number of points
    prs=max(count_data)+10 # prior took high
    nc = 6 # number of processor put high 
    fit_data = {"N":N, "d":count_data, "pes":prs, "pls":prs} # fit the data
    model_cp = stan_utility.compile_model('STAN/poisson_cp_ppc.stan')
    fit_cp = model_cp.sampling(data=fit_data, iter=R, chains=nc,seed=4838282, refresh=0) # fit the change point model to the data
    #convert fitting result to PDF curve
    bins = np.linspace(1, N, N+1)
    counts = [np.histogram(x, bins=bins)[0] for x in fit_cp["cp"][np.newaxis,:]]
    probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    counts=[c/(R*3) for c in counts]
    creds = [np.percentile([count[b] for count in counts], probs) for b in range(N)]
 
    idxs = [idx for idx in range(N)]
    xs = [bins[idx] + delta for idx in range(N) for delta in [0, bins[1]-bins[0]]]
    pad_creds = [creds[idx] for idx in idxs]
    proba=[c[4] for c in pad_creds]
    return proba



def plot_nocp (count_data,R) :
    """ Fit the no change point model on the data and plot the result. """

    N=len(count_data)#number of points
    prs=max(count_data)+10 # prior took high
    nc = 6 # number of processor put high 

    #fitting
    fit_data = {"N":N, "d":count_data, "prs":prs}
    model_nocp = stan_utility.compile_model('STAN/poisson_ppc.stan')

    # posterior predictive distribution
    Xs= fit_nocp['d_ppc']
    cumsums = np.array([np.cumsum(x) for x in Xs])
    probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    creds = [np.percentile(xt, probs) for xt in cumsums.T]
    idxs = [ idx for idx in range(Xs.shape[1]) for r in range(2) ]
    xs = [ idx + delta for idx in range(1, Xs.shape[1]+1) for delta in [-0.5,0.5]]
    pad_creds = [ creds[idx] for idx in idxs ]

    #plotting
    plt.figure()
    plt.fill_between(xs, [c[0] for c in pad_creds], [c[8] for c in pad_creds],
                    facecolor=light, color=light,label='80% predictive distribution')
    plt.fill_between(xs, [c[0] for c in pad_creds], [c[7] for c in pad_creds],
                    facecolor=light_highlight, color=light_highlight,label='60%')
    plt.fill_between(xs, [c[2] for c in pad_creds], [c[6] for c in pad_creds],
                    facecolor=mid, color=mid,label='40%')
    plt.fill_between(xs, [c[3] for c in pad_creds], [c[5] for c in pad_creds],
                    facecolor=mid_highlight, color=mid_highlight,label='20%')
    plt.plot(xs, [c[4] for c in pad_creds], color=dark,label='model',linewidth=1)
    plt.legend(bbox_to_anchor=(0.025, 1.06), loc=2, borderaxespad=0.,fontsize=12)
    plt.gca().set_xlim([0, Xs.shape[1]])
    plt.gca().set_xlabel("Data point",size=15)
    plt.xticks(fontsize=15)
    plt.gca().set_ylim([0, max([c[8] for c in creds])])
    plt.gca().set_ylabel("Cumulative number of events",size=15,color='crimson')
    plt.yticks(fontsize=12)
    plt.show()
    return



def plot_cp (count_data,R) :
    """ Fit the change point model to the data and plot the result. """

    N=len(count_data)#number of points
    prs=max(count_data)+10 # prior took high
    nc = 6 # number of processor put high 

    #fitting
    fit_data = {"N":N, "d":count_data, "pes":prs, "pls":prs}
    model_cp = stan_utility.compile_model('STAN/poisson_cp_ppc.stan')
    fit_cp = model_cp.sampling(data=fit_data, iter=R, chains=nc,seed=4838282, refresh=0)

    #probability density function
    bins = np.linspace(1, N, N+1)
    counts = [np.histogram(x, bins=bins)[0] for x in fit_cp["cp"][np.newaxis,:]]
    probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    counts=[c/(R*3) for c in counts]
    creds = [np.percentile([count[b] for count in counts], probs) for b in range(N)]
    idxs = [idx for idx in range(N)]
    xs = [bins[idx] + delta for idx in range(nbins) for delta in [0, bins[1]-bins[0]]]
    pad_creds = [creds[idx] for idx in idxs]
    proba=[c[4] for c in pad_creds]

    # posterior predictive distribution
    Xs= fit_cp['d_ppc']
    cumsums = np.array([np.cumsum(x) for x in Xs])
    probs = [10, 20, 30, 40, 50, 60, 70, 80, 90]
    creds = [np.percentile(xt, probs) for xt in cumsums.T]
    idxs = [ idx for idx in range(Xs.shape[1]) for r in range(2) ]
    xs = [ idx + delta for idx in range(1, Xs.shape[1]+1) for delta in [-0.5,0.5]]
    pad_creds = [ creds[idx] for idx in idxs ]

    #plotting
    plotline=np.cumsum(count_data)
    fig, ax1 = plt.subplots()
    ax1.fill_between(xs, [c[0] for c in pad_creds], [c[8] for c in pad_creds],
                facecolor=light, color=light,label='80% predictive distribution')
    ax1.fill_between(xs, [c[0] for c in pad_creds], [c[7] for c in pad_creds],
                facecolor=light_highlight, color=light_highlight,label='60%')
    ax1.fill_between(xs, [c[2] for c in pad_creds], [c[6] for c in pad_creds],
                facecolor=mid, color=mid,label='40%')
    ax1.fill_between(xs, [c[3] for c in pad_creds], [c[5] for c in pad_creds],
                facecolor=mid_highlight, color=mid_highlight,label='20%')
    ax1.plot(xs, [c[4] for c in pad_creds], color=dark,label='model',lw=1)
    ax1.plot(range(1, len(plotline)+1), plotline, c='k', lw=2,label='data')
    ax1.set_xlabel('Data point',size=15)
    ax1.set_ylabel('Cumulative nomber', color='crimson',size=15)
    for tl in ax1.get_yticklabels():
        tl.set_color('crimson')
    ax2 = ax1.twinx()
    ax2.plot(proba,color='navy',linewidth=1)
    ax2.set_ylim([0, 0.7])
    ax2.set_ylabel('Probability Density Function of the change point',color="navy",size=12)
    for tl in ax2.get_yticklabels():
        tl.set_color("navy")
    ax1.legend(bbox_to_anchor=(0.025, 1.06), loc=2, borderaxespad=0.,fontsize=12)
    plt.xticks(fontsize=15)
    plt.show()
    return



def lpd(count_data,R):
    """ Compute lpd values of the 2 models, keep its difference, its standard deviation and the location of the change point found. """

    N=len(count_data)#number of points
    prs=max(count_data)+10 # prior took high
    nc = 6 # number of processor put high 

    #compute lpd value of the no change model
    lpd_poisson = np.zeros (N) # list of lpd values of points
    model_poisson_loo = stan_utility.compile_model('STAN/poisson_loo.stan')
    for n in range(N): #iterate
        d_loo = count_data[n]#point removed
        d = np.array([count_data[i] for i in range(N) if i != n])#rest of data
        fit_data = {"N":N, "d":d, "prs":prs, "d_loo":d_loo, "n_loo":n}
        fit_nocp_loo = model_poisson_loo.sampling(data=fit_data,iter=R, chains=nc,seed=4838282, refresh=0)#4000 fits
        lpd_poisson[n] = logsumexp(fit_nocp_loo["log_lik_d_loo"]) - np.log(len(fit_nocp_loo["log_lik_d_loo"]))
    loonochange = np.sum(lpd_poisson) # total lpd value of the model

    # compute lpd value of the change point model
    lpd_poisson_cp = np.zeros(N) # list of lpd values of points
    model_poisson_loo_cp = stan_utility.compile_model('STAN/poisson_cp_loo.stan')
    for n in range(N):#iterate
        d_loo = count_data[n]#point removed
        d = np.array([count_data[i] for i in range(N) if i != n])#rest of data
        fit_data = {"N":N, "d":d, "pes":prs, "pls":prs, "d_loo":d_loo, "n_loo":n}
        fit_cp_loo = model_poisson_loo_cp.sampling(data=fit_data,iter=R, chains=nc,seed=4838282,refresh=0)  
        lpd_poisson_cp[n] = logsumexp(fit_cp_loo["log_lik_d_loo"]) - np.log(len(fit_cp_loo["log_lik_d_loo"]))
    loochange=np.sum(lpd_poisson_cp) # total lpd value of the model

    #compute standard deviation of the difference of lpd values
    diff = lpd_poisson_cp - lpd_poisson
    sigma = np.sqrt(N*st.pvariance(diff))

    #look for change point location
    cpdistrib = cpsearch (count_data, R)
    cp = np.argmax(cpdistrib)

    return (loochange-loonochange,sigma,cp) #lpd value difference, its standard deviation, and change point location









