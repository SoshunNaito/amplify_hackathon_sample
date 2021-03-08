import numpy as np
import random
import itertools
import queue

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

client = FixstarsClient()
client.token = "Xiccn8dKHhDoboWnaixrUEDRjvMl2vzo"
client.parameters.timeout = 1000

N = 6	# number of qubits
M = 10	# number of CX-gate layers

def calcCost(A, B): # count the number of CX-gates between 2 layers
	if(len(A) != len(B)): return -1
	_A = sorted(A)
	_B = sorted(B)

	n = len(A)

	for i in range(n):
		if(_A[i] != _B[i]): return -1

	for i in range(n-1):
		if(_A[i] == _A[i+1]): return -1
		if(_B[i] == _B[i+1]): return -1

	C = np.zeros(n)
	for i in range(n):
		for j in range(n):
			if(A[i] == B[j]):
				C[i] = j
				break

	print(C)
	ans = 0
	for i in range(n):
		for j in range(i+1,n):
			if(C[i] > C[j]):
				ans += 1
	return ans
"""
A = np.array([1,2,5,0,3,4])
B = np.array([5,3,2,4,0,1])

print(calcCost(A,B))
"""

def gen_CXgate_layers(): # generate CX-gates for each layer
	ans = []
	for m in range(M):
		n = random.randrange(2,N+1,2)
		g = random.sample(range(N), n)
		
		for i in range(len(g)):
			g[i] = chr(ord("a")+g[i])
		
		random.shuffle(g)
		ans.append(g)

	return ans

query = gen_CXgate_layers()
print(query)

RUN_CLASSIC_SOLVER = True
RUN_QUANTUM_SOLVER = True

if(RUN_CLASSIC_SOLVER):
	##########   generate all permutations   ##########
	x = []
	for v in itertools.permutations(range(N), N):
		s = ""
		for c in v:
			s += chr(ord("a")+c)
		x.append(s)
	#print(x)

	x_inv = {}
	for i in range(len(x)):
		x_inv[x[i]] = i
	#print(x_inv)


	##########   DP solution in O((N!)^2 * M)   ##########
	F = len(x)
	INF = 10000000
	dp = [[INF for f in range(F)] for m in range(M)]
	back = [[-1 for f in range(F)] for m in range(M)]

	for m in range(M):
		if(m == 0):
			for f in range(F):
				dp[m][f] = 0

				valid = True
				for g0 in range(0,len(query[m]),2):
					c0, c1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f].find(c0)-x[f].find(c1)) != 1):
						valid = False
						break

				if(valid == False):
					dp[m][f] = INF
					continue
		else:
			for f2 in range(F):
				dp[m][f2] = dp[m-1][f2]
				back[m][f2] = f2

				valid = True
				for g0 in range(0,len(query[m]),2):
					c0, c1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f2].find(c0)-x[f2].find(c1)) != 1):
						valid = False
						break

				if(valid == False):
					dp[m][f2] = INF
					continue

				# 01-BFS
				dist = [-1 for f in range(F)]
				dist[f2] = 0

				que = queue.Queue()
				que.put(f2)

				while not que.empty():
					f1 = que.get()
					s = x[f1]
					for i in range(N-1):
						c0, c1 = s[i], s[i+1]
						s = s.translate(str.maketrans({c0:c1, c1:c0}))

						f = x_inv[s]
						if(dist[f] == -1):
							dist[f] = dist[f1] + 1
							que.put(f)
							if(dp[m-1][f] + dist[f] < dp[m][f2]):
								back[m][f2] = f
								dp[m][f2] = dp[m-1][f] + dist[f]

						s = s.translate(str.maketrans({c0:c1, c1:c0}))

	cost, last = INF, -1
	for f in range(F):
		if(dp[M-1][f] < cost):
			cost, last = dp[M-1][f], f

	ans = [last] * M
	for m in range(M-1)[::-1]:
		ans[m] = back[m+1][ans[m+1]]
	ans = [x[a] for a in ans]

	print(ans)
	print("cost(classic) = " + str(cost))

if(RUN_QUANTUM_SOLVER):
	exit(0)








