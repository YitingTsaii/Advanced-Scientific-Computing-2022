import numpy as np
import random
import matplotlib.pyplot as plt
random.seed(123)

# HW2 Q1 (a)
num_traj = 10
N = 500
beta = 0.5
mu = 0.05
gamma = 0.05
tend = 150
steps = [] # to store the number of steps

for i in range(num_traj):
    S_list = []
    E_list = []
    I_list = []
    R_list = []
    t_list = []

    I = int(N*0.025)
    S = N - I
    E = 0
    R = 0
    t = 0
    S_list.append(S)
    E_list.append(E)
    I_list.append(I)
    R_list.append(R)
    t_list.append(t)
    
    while True:
        # a1, a2, a3 are propensities
        a1 = (beta/N) * S * I
        a2 = mu * E
        a3 = gamma * I
        a0 = a1 + a2 + a3
        #print("a1 =", a1, ", a2 =", a2, ", a3 = ", a3)
        if a0 == 0: # no reaction will happen any more
            S_list.append(S)
            E_list.append(E)
            I_list.append(I)
            R_list.append(R)
            t_list.append(tend)
            break
        t += random.expovariate(a0)
        if t > tend:
            t = tend
        u = random.uniform(0, a0)
        if u < a1:
            S -= 1
            E += 1
        elif u >= a1 and u < a1+a2:
            E -= 1
            I += 1
        else:
            I -= 1
            R += 1
        S_list.append(S)
        E_list.append(E)
        I_list.append(I)
        R_list.append(R)
        t_list.append(t)
        if t >= tend:
            break
        
    steps.append(len(t_list)) # record number of steps
    plt.step(t_list, S_list, where="post", color = "red")
    plt.step(t_list, E_list, where="post", color = 'blue')
    plt.step(t_list, I_list, where="post", color = 'green')
    plt.step(t_list, R_list, where="post", color = 'tab:orange')
    
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend(["S", "E", "I", "R"])

# HW2 Q1 (c)
from scipy.integrate import solve_ivp
# v = [S, E, I, R]
# S = v[0], E = v[1], I = v[2], R = v[3]
def rhs(s, v): 
    return [-(beta/N)*v[0]*v[2], (beta/N)*v[0]*v[2]-mu*v[1], mu*v[1]-gamma*v[2], gamma*v[2]]

N = 500
beta = 0.5
mu = 0.05
gamma = 0.05
tend = 150
I0 = int(N*0.025)
S0 = N - I0
E0 = 0
R0 = 0

res = solve_ivp(rhs, (0, tend), [S0, E0, I0, R0]) # (fun, time, initial_value)

t_list = res["t"]
S_list = res["y"][0]
E_list = res["y"][1]
I_list = res["y"][2]
R_list = res["y"][3]

plt.plot(t_list, S_list, color = "red")
plt.plot(t_list, E_list, color = 'blue')
plt.plot(t_list, I_list, color = 'green')
plt.plot(t_list, R_list, color = 'tab:orange')
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend(["S", "E", "I", "R"])

# HW2 Q1 (d)
num_traj = 50
N = 2000
beta = 0.5
mu = 0.05
gamma = 0.05
tend = 150

for i in range(num_traj):
    S_list = []
    E_list = []
    I_list = []
    R_list = []
    t_list = []

    I = 1
    S = N - I
    E = 0
    R = 0
    t = 0

    S_list.append(S)
    E_list.append(E)
    I_list.append(I)
    R_list.append(R)
    t_list.append(t)
    
    while True:
        a1 = (beta/N) * S * I
        a2 = mu * E
        a3 = gamma * I
        a0 = a1 + a2 + a3
        #print("a1 =", a1, ", a2 =", a2, ", a3 = ", a3)
        #print("S =", S, ", E =", E, ", I = ", I, ", R = ", R)
        if a0 == 0:
            S_list.append(S)
            E_list.append(E)
            I_list.append(I)
            R_list.append(R)
            t_list.append(tend)
            break
        t += random.expovariate(a0)
        if t >= tend:
            break
        u = random.uniform(0, a0)
        if u < a1:
            S -= 1
            E += 1
        elif u >= a1 and u < a1+a2:
            E -= 1
            I += 1
        else:
            I -= 1
            R += 1
        S_list.append(S)
        E_list.append(E)
        I_list.append(I)
        R_list.append(R)
        t_list.append(t)

    #lighter color (for comparison with ODE)
    plt.step(t_list, S_list, where="post", color = "lightcoral", linestyle = 'dotted')
    plt.step(t_list, E_list, where="post", color = 'cornflowerblue', linestyle = 'dotted')
    plt.step(t_list, I_list, where="post", color = 'mediumaquamarine', linestyle = 'dotted')
    plt.step(t_list, R_list, where="post", color = 'gold', linestyle = 'dotted')
    
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
#plt.legend(["S", "E", "I", "R"],loc='center left')

# HW2 Q1 (d) -ODE
def rhs(s, v): 
    return [-(beta/N)*v[0]*v[2], (beta/N)*v[0]*v[2]-mu*v[1], mu*v[1]-gamma*v[2], gamma*v[2]]

N = 2000
beta = 0.5
mu = 0.05
gamma = 0.05
tend = 150
I0 = 1
S0 = N - I0
E0 = 0
R0 = 0

res = solve_ivp(rhs, (0, tend), [S0, E0, I0, R0]) # (fun, time, initial_value)

t_list = res["t"]
S_list = res["y"][0]
E_list = res["y"][1]
I_list = res["y"][2]
R_list = res["y"][3]

plt.plot(t_list, S_list, color = "red", linewidth = '2')
plt.plot(t_list, E_list, color = 'blue', linewidth = '2')
plt.plot(t_list, I_list, color = 'green', linewidth = '2')
plt.plot(t_list, R_list, color = 'tab:orange', linewidth = '2')
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend(["S", "E", "I", "R"])

# HW2 Q1 (e)
num_traj = 10
tau = 10
N = 500
beta = 0.5
mu = 0.05
gamma = 0.05
tend = 150

for i in range(num_traj):
    S_list = []
    E_list = []
    I_list = []
    R_list = []
    t_list = []
    
    I = int(N*0.025)
    S = N - I
    E = 0
    R = 0
    t = 0
    
    S_list.append(S)
    E_list.append(E)
    I_list.append(I)
    R_list.append(R)
    t_list.append(t)
    
    while True:
        a1 = (beta/N) * S * I
        a2 = mu * E
        a3 = gamma * I
        #a0 = a1 + a2 + a3
        #print("a1 =", a1, ", a2 =", a2, ", a3 = ", a3)
        t += tau
        if t >= tend:
            t = tend
            #break
        k1 = np.random.poisson(a1*tau)
        k2 = np.random.poisson(a2*tau)
        k3 = np.random.poisson(a3*tau)
        # update S
        if S - k1 < 0:
            k1 = S
        S = S - k1
        # update E
        if E + k1 - k2 < 0:
            k2 = E + k1
        E = E + k1 - k2
        # update I
        if I + k2 - k3 < 0:
            k3 = I + k2
        I = I + k2 - k3
        R = R + k3
        
        S_list.append(S)
        E_list.append(E)
        I_list.append(I)
        R_list.append(R)
        t_list.append(t)
        if t >= tend:
            break
        
    steps = len(t_list) # record number of steps
    plt.step(t_list, S_list, where="post", color = "red")
    plt.step(t_list, E_list, where="post", color = 'blue')
    plt.step(t_list, I_list, where="post", color = 'green')
    plt.step(t_list, R_list, where="post", color = 'tab:orange')
    
plt.xlabel("Time (days)")
plt.ylabel("Number of individuals")
plt.legend(["S", "E", "I", "R"])

print(steps)