import numpy as np
import random
import itertools
import queue
import math

from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

client = FixstarsClient()
# client.token = "DELETED TOKEN"
client.token = "Xiccn8dKHhDoboWnaixrUEDRjvMl2vzo"
client.parameters.timeout = 10 * 1000

N = 6	# number of qubits
SAMPLE_SIZE_MAX = 100000
pairList = []
pairList_inv = {}

for i in range(N):
	for j in range(i, N-i):
		pairList_inv[i*N+j] = len(pairList)
		pairList.append(i*N+j)

for i in range(N):
	for j in range(N):
		a,b,c,d = i*N+j, j*N+i, (N-1-i)*N+(N-1-j), (N-1-j)*N+(N-1-i)
		pairList_inv[i*N+j] = pairList_inv[min(a,b,c,d)]


# print(pairList)
# print(pairList_inv)


def calcCost(A): # count the number of swap gates
	n = len(A)
	ans = 0
	for i in range(n):
		for j in range(i+1,n):
			if(A[i] > A[j]):
				ans += 1
	return ans
"""
A = np.array([1,2,5,0,3,4])
B = np.array([5,3,2,4,0,1])

print(calcCost(A,B))
"""

def getDataset(N):
	X, Y = [], []

	if(N > 15 or math.factorial(N) > SAMPLE_SIZE_MAX):
		for n in range(SAMPLE_SIZE_MAX):
			p = np.random.permutation(np.arange(N))

			x = np.zeros(len(pairList))
			for i in range(N):
				x[pairList_inv[p[i]*N + i]] += 1
			X.append(x)
			Y.append(calcCost(p))
	else:
		for p in itertools.permutations(range(N), N):
			x = np.zeros(len(pairList))
			for i in range(N):
				x[pairList_inv[p[i]*N + i]] += 1
			X.append(x)
			Y.append(calcCost(p))

	# print(X)
	# print(Y)

	return X,Y

RUN_CLASSIC_SOLVER = True
RUN_QUANTUM_SOLVER = True
PLOT_RESULTS = True

X,Y = getDataset(N)

if(RUN_CLASSIC_SOLVER):
	model = LinearRegression()
	model.fit(X,Y)

	if(PLOT_RESULTS):
		Z = model.predict(X)
		plt.scatter(Y,Z)
		plt.xlim(-1, N*(N-1)//2+1)
		plt.ylim(-1, N*(N-1)//2+1)
		plt.show()
	#print(pairList)
	print("coefficient")
	print(model.coef_)
	print("intercept")
	print(model.intercept_)

if(RUN_QUANTUM_SOLVER):
	exit(0)




