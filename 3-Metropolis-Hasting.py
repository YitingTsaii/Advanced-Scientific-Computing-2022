from scipy.stats import multivariate_normal
from scipy.stats import cauchy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# HW2 Q3
rho = 99/100
gamma = 0.5 # specify the scale parameter for cauchy
x0 = 0
y0 = 0
p = multivariate_normal(mean=[0,0], cov=[[1,rho],[rho,1]])
#sample_list = []
#sample_list.append([x0,y0])
reject_count = 0
samples = np.zeros((1000,2))
#p.pdf([1,0])

# sample from cauchy
# loc is the median, scale is the gamma
for i in range(1000):
    x1 = cauchy.rvs(loc=x0, scale=gamma, size=1)[0]
    y1 = cauchy.rvs(loc=y0, scale=gamma, size=1)[0]
    alpha = p.pdf([x1,y1])/p.pdf([x0,y0])
    
    if alpha >= 1:
        x0 = x1
        y0 = y1
    else:
        u = np.random.uniform(low=0.0, high=1.0, size=1)[0]
        if u < alpha:
            x0=x1
            y0=y1
        else:
            reject_count += 1
    #sample_list.append([x0,y0])
    samples[i][0] = x0
    samples[i][1] = y0
   
# plot samples and report rejection rates  
reject_rate = reject_count/1000
plt.scatter(samples[:,0], samples[:,1], c ="blue",s=5)
plt.xlabel("x")
plt.ylabel("y")