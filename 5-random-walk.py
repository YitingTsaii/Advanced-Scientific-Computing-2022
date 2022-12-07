import numpy as np
import random
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import math
import pandas as pd

# HW2 Q2 (a)
N = 10
M = 5
D = 2
path = np.zeros((M,N,D))
#path[M][N][D]
# e.g.
#path[0][0][1] = 1
#path[0][8][1] = 2
#path[4][9][1] = 3

# Return a random integer N such that a <= N <= b. 
# sample a dimension

for m in range(M):
    for n in range(N): # n=0 is storing the first step, the origin (0,0) is not stored
        direction = random.randint(0, D-1)
        change = random.sample(set([1, -1]), 1)[0] # either +1 or -1
        if n>0:
            for d in range(D):
                path[m][n][d] = path[m][n-1][d]
        path[m][n][direction] += change
        
# wrap this in a function
def get_path(M, N, D):
    path = np.zeros((M,N,D))
    for m in range(M):
        for n in range(N): # n=0 is storing the first step, the origin (0,0) is not stored
            direction = random.randint(0, D-1)
            change = random.sample(set([1, -1]), 1)[0] # either +1 or -1
            if n>0:
                for d in range(D):
                    path[m][n][d] = path[m][n-1][d]
            path[m][n][direction] += change
    return path

# to use
get_path(M, N, D)

# HW2 Q2 (b)
# get Euclidean distance in D dimension
def ecucliden_dis(my_array):  # distance from the origin
    return np.sqrt(np.sum(np.square(my_array)))
#ecucliden_dis(path[0][9][:])

N = 1000
M = 100
D = 2

my_path = get_path(M, N, D)

n_list = []
d_bar_list = []
for n in range(N):
    dist_M_samples = []
    for m in range(M):
        dist_M_samples.append(ecucliden_dis(my_path[m][n][:]))
    n_list.append(n+1)
    d_bar_list.append(sum(dist_M_samples)/len(dist_M_samples))
    
# plot d_bar against n
plt.plot(n_list, d_bar_list)
plt.xlabel("n")
plt.ylabel("Average Euclidean distance")

# linear regression
# log(d_bar) = log(b) + a*log(n)
X = np.array(np.log(n_list)).reshape(-1, 1) 
y = np.array(np.log(d_bar_list))
reg = LinearRegression().fit(X, y)
reg.score(X, y)
a = reg.coef_[0]
b = math.exp(reg.intercept_)

#reg.predict(X)
print("a =", a)
print("b =", b)

pred_list = [b*(n**a) for n in n_list]

# log scale plot
plt.plot(n_list, d_bar_list, color = 'blue')
plt.plot(n_list, pred_list, color = 'red', linestyle = "dashed")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("log(n)")
plt.ylabel("log(Euclidean distance)")
plt.legend(["simulated", "fitted"])

# normal scale plot
plt.plot(n_list, d_bar_list, color = 'blue')
plt.plot(n_list, pred_list, color = 'red', linestyle = "dashed")
plt.xlabel("n")
plt.ylabel("Euclidean distance")
plt.legend(["simulated", "fitted"])

# HW2 Q2 (c)
# need to run (b) beforehand
R = 10
passing_time = np.zeros((M,R))
#passing_time[m][r] = 1

for m in range(M):
    for r in range(1, R+1): #r=1,2,3,...,10
        for n in range(N):
            dist = ecucliden_dis(my_path[m][n][:])
            if dist >= r:
                passing_time[m][r-1] = n+1
                break

passing_time_df = pd.DataFrame(passing_time)
time_mean = passing_time_df.mean()
time_std = passing_time_df.std()
r_list = [i for i in range(1,R+1)]

# linear regression
# log(d_bar) = log(b) + a*log(n)
X = np.array(np.log(r_list)).reshape(-1, 1) 
y = np.array(np.log(time_mean))
reg = LinearRegression().fit(X, y)
reg.score(X, y)
a = reg.coef_[0]
b = math.exp(reg.intercept_)

#reg.predict(X)
print("a =", a)
print("b =", b)

pred_list = [b*(r**a) for r in r_list]

# plot mean + error bar
plt.plot(r_list, time_mean, color = 'blue')
plt.errorbar(x = r_list, y = time_mean, yerr = time_std, fmt ='o')
plt.plot(r_list, pred_list, color = 'red', linestyle = "dotted")
plt.xlabel("R")
plt.ylabel("First passage time")
plt.legend(["simulated", "fitted"], loc='upper left')

plt.plot(r_list, time_mean, color = 'blue')
plt.plot(r_list, pred_list, color = 'red', linestyle = "dashed")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("log(R)")
plt.ylabel("log(First passage time)")
plt.legend(["simulated", "fitted"])