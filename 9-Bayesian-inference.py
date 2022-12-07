import pandas as pd
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal
from scipy.stats import norm # for 1D normal density
from scipy.optimize import fsolve
import matplotlib.lines as mlines

path = 'data.csv'
data = pd.read_csv(path)

# (a)
x = [i for i in data['x']]
y = [i for i in data['y']]
xy = [i*j for i,j in zip(x,y)]
x2 = [i*i for i in x]
a_MLE = sum(xy)/sum(x2)
ax_y_2 = [(a_MLE*i-j)**2 for i,j in zip(x,y)]
sigma_MLE = (sum(ax_y_2)/len(x))**0.5
print('a_MLE =', a_MLE)
print('sigma_MLE =', sigma_MLE)
# check the second order derivative is negative definite
m = [(a_MLE*i-j)*i for i,j in zip(x,y)]
l = [(a_MLE*i-j)**2 for i,j in zip(x,y)]
N = len(x)
upper_left = -sum(x2) / (sigma_MLE**2)
upper_right = 2*sum(m) / (sigma_MLE**3) 
lower_left = upper_right
lower_right = N/(sigma_MLE**2) - 3*sum(l)/(sigma_MLE**4) 
sec_der = [[upper_left, upper_right],
           [lower_left, lower_right]]
print()
print("Check the second order derivative is negative definite")
print('a =', upper_left)
print('a < 0? => ', upper_left <0 )
print('determinant =', upper_left*lower_right-upper_right*lower_left)
print('determinant > 0? => ', upper_left*lower_right-upper_right*lower_left > 0)

# (b): the Laplace approximate posterior is N(theta_MAP, cov_mat)
# get the theta_MAP
a_MAP = sum(xy)/sum(x2)
ax_y_2 = [(a_MAP*i-j)**2 for i,j in zip(x,y)]
sigma_MAP = (sum(ax_y_2)/len(x))**0.5
theta_MAP = [a_MAP, sigma_MAP]

# get the Sigma covariance matrix (for the approximated 2D Gaussian)
# second derivative of L
m = [(a_MAP*i-j)*i for i,j in zip(x,y)]
l = [(a_MAP*i-j)**2 for i,j in zip(x,y)]
N = len(x)
upper_left = -sum(x2) / (sigma_MAP**2)
upper_right = 2*sum(m) / (sigma_MAP**3) 
lower_left = upper_right
lower_right = N/(sigma_MAP**2) - 3*sum(l)/(sigma_MAP**4) 
neg_sec_der = [[-upper_left, -upper_right],
               [-lower_left, -lower_right]]
cov_mat = np.linalg.inv(neg_sec_der)

# (c)
def posterior(a, sigma):
    ans = 1
    for i in range(N):
        current = (1/(2*np.pi*(sigma**2))**0.5) * np.exp(-(a*x[i]-y[i])**2/(2*sigma**2))
        ans = ans*current
    return ans

def log_posterior(a, sigma):
    k = [(a*i-j)**2 for i,j in zip(x,y)]
    return -(N/2)*np.log(2*np.pi) - N*np.log(sigma) - sum(k)/(2*(sigma**2)) + const

# sample with Metropolis-Hasting
# change this so that the rejection rate is around 0.5
scale = 0.1 # specify the scale parameter for proposed distribution
x0 = a_MAP # set this as the initial value
y0 = sigma_MAP # set this as the initial value
#p = multivariate_normal(mean=theta_MAP, cov=cov_mat) # Goal: sample from this
#sample_list = []
#sample_list.append([x0,y0])
reject_count = 0
samples = np.zeros((1000,2))
#p.pdf([1,0])

# sample from normal
# loc is the median, scale is the gamma
for i in range(1000):
    x1 = random.normal(loc=x0, scale=scale, size=1)[0]
    y1 = random.normal(loc=y0, scale=scale, size=1)[0]
    #alpha = p.pdf([x1,y1])/p.pdf([x0,y0])
    alpha = posterior(x1,y1)/posterior(x0,y0)
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
# plot the samples by themselves
reject_rate = reject_count/1000
print("rejection rate =", reject_rate)
plt.scatter(samples[:,0], samples[:,1], c ="blue",s=5)
plt.xlabel("a")
plt.ylabel("sigma")


# plot the contours
const = 0
x = [i for i in data['x']]
y = [i for i in data['y']]
N = len(x)

x_pts = np.linspace(0, 3, 100) # a
y_pts = np.linspace(0.01, 5, 100) # sigma
X, Y = np.meshgrid(x_pts, y_pts)
fig, ax = plt.subplots()
# log posterior
Z = log_posterior(X, Y)
plt.contour(X, Y, Z, levels = np.linspace(-50,0,20), colors='black',linestyles="-")
plt.contour(X, Y, Z, levels = [-200, -150,-100, -80, -70, -60, -55], colors='black',label='Log Posterior',linestyles="-")
pos = np.dstack((X,Y))
rv = multivariate_normal(theta_MAP, cov_mat)
# log Laplace
log_la = np.log(rv.pdf(pos))
plt.contour(X, Y, log_la, levels = [-200, -150,-100, -80, -70, -60, -55,-20,-10,-5,-2,-0.01],colors='green', label='Laplace approx',linestyles="-")
# add the samples
plt.scatter(samples[:,0], samples[:,1], c ="blue",s=1, label='samples', alpha=0.3)
# add the MAP
plt.scatter(a_MAP, sigma_MAP, label='MAP', c='orange')
# add the legend and labels
plt.xlabel("a")
plt.ylabel("sigma")
#plt.legend()

my_legend1 = mlines.Line2D([], [], color='black',markersize=15, label='log-posterior')
my_legend2 = mlines.Line2D([], [], color='green',markersize=15, label='log-Laplace approximation')
my_legend3 = mlines.Line2D([], [], color='white',marker='o',markerfacecolor='orange',markersize=8, label='MAP')
ax.legend(handles=[my_legend1, my_legend2, my_legend3])
my_legend4 = mlines.Line2D([], [], color='white',marker='o',markerfacecolor='blue',markersize=6, label='samples')
ax.legend(handles=[my_legend1, my_legend2, my_legend3, my_legend4])
plt.show()


# (d)
# "samples" are the 1000 samples from the posterior using  MH
samples.mean(axis=0)  # see the samples mean (for a and sigma)
a_mean = samples.mean(axis=0)[0]

def normal_pdf(y, mean, sigma):
    # a functino of y
    # this is exactly the same as norm.pdf(y, loc, scale)
    return (2 * np.pi * sigma**2)**(-0.5) * np.exp(-(y-mean)**2/(2*sigma**2))

def cdf_P(y, x, samples):
    # a function of y
    K = len(samples)
    cdf_sum = 0
    for i in range(K):
        a_k = samples[i,0]
        sigma_k = samples[i,1]
        cdf_sum += norm.cdf(y, loc=a_k*x, scale=sigma_k)
    return cdf_sum/K

# inverse the CDF, solve by fsolve
def CDF_inv(q,x,samples):
    def solve0(y):
        return cdf_P(y, x=x, samples=samples) - q
    y_ans = fsolve(solve0, a_mean*x) # initial value set as the regression mean
    return y_ans
#try
#CDF_inv(q=0.95, x=-5, samples=samples)[0]
    
df = pd.DataFrame(columns=['x','lower90','lower50','median','upper50', 'upper90'])
xs = np.linspace(-5,5,200)
#len(xs) #200
for x in xs:
    median = CDF_inv(q=0.5, x=x, samples=samples)[0]
    lower90 = CDF_inv(q=0.05, x=x, samples=samples)[0]
    upper90 = CDF_inv(q=0.95, x=x, samples=samples)[0]
    lower50 = CDF_inv(q=0.25, x=x, samples=samples)[0]
    upper50 = CDF_inv(q=0.75, x=x, samples=samples)[0]
    df1 = pd.DataFrame([[x,lower90, lower50, median, upper50, upper90]], columns=['x','lower90', 'lower50', 'median', 'upper50', 'upper90'])
    df = pd.concat([df, df1])
    print(x)

path = 'Q1d_df.csv'
df.to_csv(path)
df = pd.read_csv(path, index_col = 0)

# plot
fig,ax = plt.subplots()
ax.plot(df['x'], df['median'], '-', color='blue', label = 'median')
ax.set_title('Prediction with Confidence Interval')
ax.set_ylabel('y')
ax.set_xlabel('x')

# plot our confidence band
ax.fill_between(df['x'], df['lower90'], df['upper90'], alpha=0.2, color='tab:orange', label = '90% CI')
ax.fill_between(df['x'], df['lower50'], df['upper50'], alpha=0.2, color='blue', label = '50% CI')
x = [i for i in data['x']]
y = [i for i in data['y']]
ax.scatter(x, y, label = 'original data')
ax.legend(loc='lower right')