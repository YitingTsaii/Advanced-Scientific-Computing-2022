# Q1 (a)
import korali
import numpy as np
import math

def Ackley(x, y):
    return -20*np.exp(-0.2 * (0.5*( x**2 + y**2 ))**0.5 )\
        - np.exp( 0.5*( np.cos(2*math.pi*x)+np.cos(2*math.pi*y) ) ) + 20 + np.exp(1)

if __name__ == '__main__':

    # The optimization problem is described in a korali Experiment object
    e = korali.Experiment()

    # Korali requires a specific interface for the function to maximize
    def objective_function(ksample):
        x, y = ksample["Parameters"]
        # note the minus sign because we want to minimize the function
        ksample["F(x)"] = -Ackley(x, y)

    e["Problem"]["Type"] = "Optimization"
    e["Problem"]["Objective Function"] = objective_function

    # Defining the problem's variables.
    e["Variables"][0]["Name"] = "x"
    e["Variables"][0]["Initial Value"] = 5.0
    e["Variables"][0]["Initial Standard Deviation"] = 5.0

    e["Variables"][1]["Name"] = "y"
    e["Variables"][1]["Initial Value"] = 5.0
    e["Variables"][1]["Initial Standard Deviation"] = 5.0

    # Configuring CMA-ES parameters
    e["Solver"]["Type"] = "Optimizer/CMAES"
    e["Solver"]["Population Size"] = 16 # lambda
    e["Solver"]["Mu Value"] = 2 # mu     
    e["Solver"]["Termination Criteria"]["Min Value Difference Threshold"] = 1e-32
    e["Solver"]["Termination Criteria"]["Max Generations"] = 100
    
    # Experiments need a Korali Engine object to be executed
    k = korali.Engine()

    # Run the optimization
    k.run(e)

# Q1 (b)
import numpy as np
from scipy.optimize import minimize
import math
import matplotlib.pyplot as plt

# x = x[0], y = x[1]
def Ackley(x):
    return -20*np.exp(-0.2 * (0.5*( x[0]**2 + x[1]**2 ))**0.5 )\
        - np.exp( 0.5*( np.cos(2*math.pi*x[0])+np.cos(2*math.pi*x[1]) ) ) + 20 + np.exp(1)

def save(loc):
    global P
    P.append(loc)

# the jacobian function
def jac(vec):
    x, y = vec
    F = [0, 0]
    a = 0.5*( x**2 + y**2 )
    b = 0.5*( np.cos(2*math.pi*x)+np.cos(2*math.pi*y) )
    F[0] = 2 * x * np.exp(-0.2*a**0.5) * a**(-0.5) + np.exp(b) * np.sin(2*math.pi*x) * math.pi
    F[1] = 2 * y * np.exp(-0.2*a**0.5) * a**(-0.5) + np.exp(b) * np.sin(2*math.pi*y) * math.pi
    return F

# minimize using python package   
x0 = np.array([5.0, 5.0])
P = [ x0 ]
minimize(Ackley, x0, method='CG', tol=1e-100, callback=save, jac=jac, options={'maxiter':100} )

i = 0
for x in P:
    print("Iteration =", i, ": (x,y) =", x, ", F(x,y) =", Ackley(x))
    i += 1
  
# plot the convergence plot (x,y)
xs = []
ys = []
for i in range(len(P)):
    xs.append(P[i][0])
    ys.append(P[i][1])
index = range(len(xs))
plt.plot(index, xs, linestyle='-', marker='o', label='x')
plt.plot(index, ys, linestyle='-', marker='o', label='y')
plt.xticks(np.arange(0, 4, 1))
plt.legend()
plt.show()

# plot the convergence plot (F)
Fs = []
for x in P:
    #print("(x,y) =", x, "F(x,y) =", Ackley(x))
    Fs.append(Ackley(x))

index = range(len(Fs))
plt.plot(index, Fs, linestyle='-', marker='o', label='Ackley(x,y)')
plt.xticks(np.arange(0, 4, 1))
plt.legend()
plt.show()
  

# vanilla gradient descent
x0 = [5, 5]
x1 = [0, 0]
gamma = 0.001 # step size
x_path = [5]
y_path = [5]
ftn_path = [Ackley(x0)]
iter = 0
while True:
    iter += 1
    grad = jac(x0)
    x1[0] = x0[0] - gamma*grad[0]  
    x1[1] = x0[1] - gamma*grad[1]  
    x_path.append(x1[0])
    y_path.append(x1[1])
    ftn_path.append(Ackley(x1))
    #loss = np.sqrt((x1[0]-x0[0])**2 + (x1[1]-x0[1])**2)
    #loss = abs((x1[0]-x0[0])) + abs((x1[1]-x0[1]))
    loss= abs(Ackley(x1)-Ackley(x0))
    #print(loss)
    if loss < 1e-100 or iter>100:
        break
    x0 = x1

count = 0
for i in range(len(x_path)):
    print("Iteration =", count, ": (x,y) =(", x_path[i], ",", y_path[i], ")," , "F(x,y) =", ftn_path[i])
    count += 1

# plot the convergence plot (x,y)
index = range(len(x_path))
plt.plot(index, x_path, linestyle='-', marker='o', label='x')
plt.plot(index, y_path, linestyle='-', marker='o', label='y')
plt.xticks(np.arange(0, 4, 1))
plt.legend()
plt.show()

# plot the convergence plot (F)
Fs = []
for x in P:
    #print("(x,y) =", x, "F(x,y) =", Ackley(x))
    Fs.append(Ackley(x))

index = range(len(ftn_path))
plt.plot(index, ftn_path, linestyle='-', marker='o', label='Ackley(x,y)')
plt.xticks(np.arange(0, 4, 1))
plt.legend()
plt.show()
