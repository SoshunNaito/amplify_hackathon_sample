import numpy as np
import random
import itertools
import math

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt

N = 6	# number of qubits
SAMPLE_SIZE_MAX = 500000
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

def func(i,j):
	return (i+j)/2 - i*j/(N-1)

def approxTest(A):
	ans = 0
	for i in range(N):
		ans += func(i, A[i])
	return ans


def getDataset(N):
	X, Y, W = [], [], []

	if(N > 15 or math.factorial(N) > SAMPLE_SIZE_MAX):
		for n in range(SAMPLE_SIZE_MAX):
			p = np.random.permutation(np.arange(N))

			x = np.zeros(len(pairList))
			for i in range(N):
				x[pairList_inv[p[i]*N + i]] += 1
			X.append(x)
			Y.append(calcCost(p))
			W.append(approxTest(p))
	else:
		for p in itertools.permutations(range(N), N):
			x = np.zeros(len(pairList))
			for i in range(N):
				x[pairList_inv[p[i]*N + i]] += 1
			X.append(x)
			Y.append(calcCost(p))
			W.append(approxTest(p))

	# print(X)
	# print(Y)

	return X,Y,W

PLOT_RAW_RESULTS = False
PLOT_RESULTS = True

X,Y,W = getDataset(N)

model_XY = LinearRegression()
model_XY.fit(X,Y)

# print(pairList)
# print("coefficient")
# print(model_XY.coef_)
# print("intercept")
# print(model_XY.intercept_)

x = np.arange(-1, N*(N-1)//2+1)

W = np.array(W).reshape(-1, 1)
Z = model_XY.predict(X)

model_WZ = LinearRegression()
model_WZ.fit(W,Z)

if(PLOT_RAW_RESULTS):
	plt.scatter(Y,Z)
	plt.xlim(-1, N*(N-1)//2+1)
	plt.ylim(-1, N*(N-1)//2+1)
	plt.xlabel("ground truth")
	plt.ylabel("sklearn results")
	plt.plot(x, x, color = "red")
	plt.show()

	plt.scatter(Y,W)
	plt.xlim(-1, N*(N-1)//2+1)
	plt.ylim(-1, N*(N-1)//2+1)
	plt.xlabel("ground truth")
	plt.ylabel("theory results")
	plt.plot(x, x, color = "red")
	plt.show()

if(PLOT_RESULTS):
	y = x * model_WZ.coef_[0] + model_WZ.intercept_

	plt.scatter(W,Z)
	plt.xlim(-1, N*(N-1)//2+1)
	plt.ylim(-1, N*(N-1)//2+1)
	plt.xlabel("theory results")
	plt.ylabel("sklearn results")
	plt.title("N = "+str(N))
	plt.plot(x,y,
		label="y = "+str(model_WZ.coef_[0])+"x"+str(model_WZ.intercept_).replace("-"," - ")+"\n"+
		"R^2 = "+str(r2_score(model_WZ.predict(W),Z)))
	plt.legend()
	plt.show()




