# Q3 - gradient based method
import random
import scipy.optimize


def fun(x):
	E = 0
	m = len(x)
	for i in range(0, m, 3):
		for j in range(i + 3, m, 3):
			dx = x[i] - x[j]
			dy = x[i + 1] - x[j + 1]
			dz = x[i + 2] - x[j + 2]
			r2 = dx**2 + dy**2 + dz**2
			r2inv = 1 / r2
			r6inv = r2inv * r2inv * r2inv
			r12inv = r6inv * r6inv
			E += r12inv - r6inv
	return 4 * E


def jac(x):
	m = len(x)
	F = [0] * m
	for i in range(0, m, 3):
		for j in range(i + 3, m, 3):
			dx = x[i] - x[j]
			dy = x[i + 1] - x[j + 1]
			dz = x[i + 2] - x[j + 2]
			r2 = dx**2 + dy**2 + dz**2
			r2inv = 1 / r2
			r4inv = r2inv * r2inv
			r8inv = r4inv * r4inv
			r14inv = r8inv * r4inv * r2inv
			f = 24 * r8inv - 48 * r14inv
			fx = dx * f
			fy = dy * f
			fz = dz * f
			F[i] += fx
			F[i + 1] += fy
			F[i + 2] += fz
			F[j] -= fx
			F[j + 1] -= fy
			F[j + 2] -= fz
	return F

random.seed(12345)
fun_value_list = []
x_list = []
n = 32 # change here
for t in range(10):
	x0 = [random.uniform(-1, 1) for i in range(3 * n)]
	res = scipy.optimize.minimize(fun, x0, method='CG', jac=jac)
	print(t, res.message, res.fun)
	print(res.x)
	fun_value_list.append(res.fun)
	x_list.append(res.x)


min_idx = fun_value_list.index(min(fun_value_list))
min_energy = fun_value_list[min_idx]
min_position = x_list[min_idx]

print("N =", n)
print("E =", min_energy)
print("Positions (xi, yi, zi) =")
for i in range(0, 3*n, 3):
    print("(", round(min_position[i],4), ", ",round(min_position[i+1],4), ", ", round(min_position[i+2],4), ")", sep='')

# Q3 - stochastic method
import math
import numpy as np
import random
import scipy.linalg

def fun(x):
	E = 0
	m = len(x)
	for i in range(0, m, 3):
		for j in range(i + 3, m, 3):
			dx = x[i] - x[j]
			dy = x[i + 1] - x[j + 1]
			dz = x[i + 2] - x[j + 2]
			r2 = dx**2 + dy**2 + dz**2
			r2inv = 1 / r2
			r6inv = r2inv * r2inv * r2inv
			r12inv = r6inv * r6inv
			E += r12inv - r6inv
	return 4 * E

def cmaes(f, x, sigma, g_max, Trace=False):

	def cumulation(c, A, B):
		alpha = 1 - c
		beta = math.sqrt(c * (2 - c) * mueff)
		return [alpha * a + beta * b for a, b in zip(A, B)]

	def wsum(A):
		return [
		    math.fsum(w * a[i] for w, a in zip(weights, A)) for i in range(N)
		]

	xmean, N = x[:], len(x)
	lambd = 4 + int(3 * math.log(N))
	mu = lambd // 2
	weights = [math.log((lambd + 1) / 2) - math.log(i + 1) for i in range(mu)]
	weights = [e / math.fsum(weights) for e in weights]
	mueff = 1 / math.fsum(e**2 for e in weights)
	cc = (4 + mueff / N) / (N + 4 + 2 * mueff / N)
	cs = (mueff + 2) / (N + mueff + 5)
	c1 = 2 / ((N + 1.3)**2 + mueff)
	cmu = min(1 - c1, 2 * (mueff - 2 + 1 / mueff) / ((N + 2)**2 + mueff))
	damps = 1 + 2 * max(0, math.sqrt((mueff - 1) / (N + 1)) - 1) + cs
	chiN = math.sqrt(2) * math.gamma((N + 1) / 2) / math.gamma(N / 2)
	ps, pc, C = [0] * N, [0] * N, np.identity(N)
	trace = []
	for gen in range(1, g_max + 1):
		sqrtC = np.real(scipy.linalg.sqrtm(C))
		x0 = [[random.gauss(0, 1) for d in range(N)] for i in range(lambd)]
		x1 = [sqrtC @ e for e in x0]
		xs = [xmean + sigma * e for e in x1]
		ys = [f(e) for e in xs]
		ys, x0, x1, xs = zip(*sorted(zip(ys, x0, x1, xs)))
		xmean = wsum(xs)
		ps = cumulation(cs, ps, wsum(x0))
		pssq = math.fsum(e**2 for e in ps)
		sigma *= math.exp(cs / damps * (math.sqrt(pssq) / chiN - 1))
		Cmu = sum(w * np.outer(d, d) for w, d in zip(weights, x1))
		if (N + 1) * pssq < 2 * N * (N + 3) * (1 - (1 - cs)**(2 * gen)):
			pc = cumulation(cc, pc, wsum(x1))
			C1 = np.outer(pc, pc)
			C = (1 - c1 - cmu) * C + c1 * C1 + cmu * Cmu
		else:
			pc = [(1 - cc) * e for e in pc]
			C1 = np.outer(pc, pc)
			C = (1 - c1 - cmu) * C + c1 * (C1 + cc * (2 - cc) * C) + cmu * Cmu
		if Trace:
			trace.append(
			    (gen * lambd, ys[0], xs[0], sigma, C, ps, pc, Cmu, C1, xmean))
	return trace if Trace else xmean

    
## use fun
random.seed(12345)
n=8 # change here
trace = cmaes(fun, [random.uniform(-1, 1) for i in range(3*n)], 0.1, 5000, Trace=True)
last_G = trace[-1]  # get the last iteration
neval, min_energy, min_position, *rest = last_G
print("N =", n)
print("E =", min_energy)
print("Positions (xi, yi, zi) =")
for i in range(0, 3*n, 3):
    print("(", round(min_position[i],4), ", ",round(min_position[i+1],4), ", ", round(min_position[i+2],4), ")", sep='')

