# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 17:20:26 2018

@author: katsuhisa

Reproducing similar figures to Sanborn & Charter, Trends in Cognitive Science 2016.
As a Bayesian sampling, a simple Metropolis-Hastings algorithm is implemented.
"""

# libraries ------------------------------
import numpy as np
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

#%%
# functions -----------------------------
def target_params(disttype):
    """preset parameters for target distributions"""
    if disttype == 'round':
        me = np.array([2,2])
        cov = np.array([[1,0],[0,1]])
    elif disttype == 'correlated':
        me = np.array([2,2])
        cov = np.array([[1,0.9],[0.9,1]])
    elif disttype == 'close_bimodal':
        me = np.array([[0,0],[4,4]])
        cov = np.array([[1,0],[0,1]])
    elif disttype == 'bimodal':
        me = np.array([[0,0],[10,10]])
        cov = np.array([[1,0],[0,1]])
    return me, cov    
        
def target_rvs(N, disttype):
    """random samples from target distribution (bivariate Gaussian)"""
    me, cov = target_params(disttype)
    if disttype in ['round', 'correlated']:
        p = multivariate_normal(mean=me, cov=cov).rvs(size=N, random_state=1220).T
        return p[0], p[1]
    elif disttype in ['bimodal', 'close_bimodal']:
        p1 = multivariate_normal(mean=me[0], cov=cov).rvs(size=round(N/2), random_state=1220).T 
        p2 = multivariate_normal(mean=me[1], cov=cov).rvs(size=N-round(N/2), random_state=1220).T 
        return np.hstack((p1[0],p2[0])), np.hstack((p1[1],p2[1]))

def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)"""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins, 
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)

def target_pdf(p, disttype):
    """target distribution (bivariate Gaussian)"""
    me, cov = target_params(disttype)
    if disttype == 'round' or disttype == 'correlated':
        prob = multivariate_normal.pdf(p, mean=me, cov=cov)
    elif disttype == 'bimodal' or disttype == 'close_bimodal':
        prob0 = multivariate_normal.pdf(p, mean=me[0], cov=cov)
        prob1 = multivariate_normal.pdf(p, mean=me[1], cov=cov)
        prob = max([prob0, prob1])        
        
    return prob

def proposal_pdf(p):
    """proposal distribution"""
    return (p[0] + np.random.normal(0, 1), p[1] + np.random.normal(0, 1))

def mh(N, disttype):
    """metropolis hastings sampling (no 'thin' and 'burn_in' for simplicity)"""
    xs = np.array([])
    ys = np.array([])
    pos_now = (0,0)
    accept = 0
    for i in range(N):
        pos_cand = proposal_pdf(pos_now)
        prob_stay = target_pdf(pos_now, disttype)
        prob_move = target_pdf(pos_cand, disttype)
        if prob_move / prob_stay > np.random.uniform(0,1,1):
            pos_now = pos_cand
            xs = np.append(xs, pos_now[0])
            ys = np.append(ys, pos_now[1])
            accept += 1
    return xs, ys, accept/N

#%%
# visualization similar to the Fig 1b in Sanborn & Charter 2016 -----------
# preset parameters
disttype = ['round', 'correlated', 'close_bimodal', 'bimodal']
axisrange = [[-2, 6], [-2, 6], [-4, 9], [-4, 14]]
plt.close('all')

# probability distributions
fig1, axes1 = plt.subplots(1, 4)
for i in range(4):
    x, y = target_rvs(5000, disttype[i])
    xx, yy, zz = kde2D(x, y, 1)
    axes1[i].pcolormesh(xx, yy, zz)
    axes1[i].set_title(disttype[i])
    axes1[i].set_xlabel('x')
    axes1[i].set_ylabel('y')
    axes1[i].get_xaxis().set_ticks([])
    axes1[i].get_yaxis().set_ticks([])
    axes1[i].get_xaxis().set_ticklabels([])
    axes1[i].get_yaxis().set_ticklabels([])

# iterations & MH samples
fig2, axes2 = plt.subplots(2, 4)
fig3, axes3 = plt.subplots(1, 4)
for i in range(4):
    # sampling
    xs, ys, accept_ratio = mh(1000, disttype[i])
    
    # iteration
    axes2[0,i].plot(np.arange(len(xs)), xs, '-ko', alpha=0.2)
    axes2[0,i].set_title(disttype[i])
    axes2[1,i].plot(np.arange(len(ys)), ys, '-ko', alpha=0.2)
    axes2[1,i].set_xlabel('iteration')
    axes2[0,i].set_ylabel('x')
    axes2[1,i].set_ylabel('y')
    axes2[0,i].get_xaxis().set_ticks([])
    axes2[0,i].get_yaxis().set_ticks([])
    axes2[0,i].get_xaxis().set_ticklabels([])
    axes2[0,i].get_yaxis().set_ticklabels([])
    axes2[0,i].set_ylim(axisrange[i][0], axisrange[i][1])
    axes2[1,i].get_xaxis().set_ticks([])
    axes2[1,i].get_yaxis().set_ticks([])
    axes2[1,i].get_xaxis().set_ticklabels([])
    axes2[1,i].get_yaxis().set_ticklabels([])
    axes2[1,i].set_ylim(axisrange[i][0], axisrange[i][1])
    
    # sampling distributions
    axes3[i].scatter(xs, ys, s=20, c='k', alpha=0.2)
    axes3[i].set_title(disttype[i])
    axes3[i].set_xlabel('x')
    axes3[i].set_ylabel('y')
    axes3[i].set_xlim(axisrange[i][0], axisrange[i][1])
    axes3[i].set_ylim(axisrange[i][0], axisrange[i][1])
    axes3[i].get_xaxis().set_ticks([])
    axes3[i].get_yaxis().set_ticks([])
    axes3[i].get_xaxis().set_ticklabels([])
    axes3[i].get_yaxis().set_ticklabels([])
    
fig2.tight_layout()
fig3.tight_layout()