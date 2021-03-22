import numpy as np
import random
import itertools
import queue

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

import time

client = FixstarsClient()
client.token = "DELETED TOKEN"
client.parameters.timeout = 1 * 1000

constraintWeight = 100

def calcCost(A, B): # count the number of swap gates between 2 layers
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

	# print(C)
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

def gen_CXgate_layers(N,M): # generate CX-gates for each layer
	ans = []
	for m in range(M):
		n = random.randrange(2,N+1,2)
		g = random.sample(range(N), n)
		
		for i in range(len(g)):
			g[i] = chr(ord("a")+g[i])
		
		random.shuffle(g)
		ans.append(g)

	return ans

def classic_solver(N,M,query): # solve with a classic computer
	t0 = time.time()
	
	##########   generate all permutations   ##########
	x = []
	for p in itertools.permutations(range(N), N):
		s = ""
		for c in p:
			s += chr(ord("a")+c)
		x.append(s)
	# print(x)

	x_inv = {}
	for i in range(len(x)):
		x_inv[x[i]] = i
	# print(x_inv)


	##########   DP solution in O((N!)^2 * M)   ##########
	F = len(x)
	INF = 10000000
	valid = [[] for m in range(M)]
	dp = [[INF for f in range(F)] for m in range(M)]
	back = [[-1 for f in range(F)] for m in range(M)]

	dist = [[-1 for f1 in range(F)] for f2 in range(F)]

	for f1 in range(F):
		dist[f1][f1] = 0
		que = queue.Queue()
		que.put(f1)

		while not que.empty():
			f2 = que.get()
			s = x[f2]
			for i in range(N-1):
				v0, v1 = s[i], s[i+1]
				s = s.translate(str.maketrans({v0:v1, v1:v0}))

				f = x_inv[s]
				if(dist[f1][f] == -1):
					dist[f1][f] = dist[f1][f2] + 1
					que.put(f)

				s = s.translate(str.maketrans({v0:v1, v1:v0}))

	for m in range(M):
		if(m == 0):
			for f in range(F):
				dp[m][f] = 0

				flag = True
				for g0 in range(0,len(query[m]),2):
					v0, v1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f].find(v0)-x[f].find(v1)) != 1):
						flag = False
						break

				if(flag == False):
					dp[m][f] = INF
					continue
				else:
					valid[m].append(f)
		else:
			for f in range(F):
				dp[m][f] = dp[m-1][f]
				back[m][f] = f

				flag = True
				for g0 in range(0,len(query[m]),2):
					v0, v1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f].find(v0)-x[f].find(v1)) != 1):
						flag = False
						break

				if(flag == False):
					dp[m][f] = INF
					continue
				else:
					valid[m].append(f)

				for f0 in valid[m-1]:
					if(dp[m-1][f0] + dist[f0][f] < dp[m][f]):
						dp[m][f] = dp[m-1][f0] + dist[f0][f]
						back[m][f] = f0

	cost, last = INF, -1
	for f in range(F):
		if(dp[M-1][f] < cost):
			cost, last = dp[M-1][f], f

	ans = [last] * M
	for m in range(M-1)[::-1]:
		ans[m] = back[m+1][ans[m+1]]
	ans = [x[a] for a in ans]

	# print("cost(classic) = " + str(cost))
	# print(ans)

	t1 = time.time()

	return (t1-t0), cost, ans

def quantum_solver(N,M,query): # solve with Amplify
	t0 = time.time()

	q = gen_symbols(BinaryPoly, M, N, N)	# represent the solution

	##########   constraints   ##########

	# each layer doesn't have 2+ same values
	one_hot_constraints_layer = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q[m][n][v] for n in range(N)), 1)
		for m in range(M)
		for v in range(N)
	]

	# each qubit doesn't have 2+ values
	one_hot_constraints_num = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q[m][n][v] for v in range(N)), 1)
		for m in range(M)
		for n in range(N)
	]

	# Every CX-gate must be applied for 2 adjacent qubits
	CXgate_constraints = []
	for m in range(M):
		for g0 in range(0,len(query[m]),2):
			v0, v1 = ord(query[m][g0])-ord("a"), ord(query[m][g0 + 1])-ord("a")

			# v0 and v1 must be adjacent each other
			for i in range(N):
				for j in range(i+2, N):
					CXgate_constraints.append(
						penalty(q[m][i][v0] * q[m][j][v1])
					)
					CXgate_constraints.append(
						penalty(q[m][i][v1] * q[m][j][v0])
					)

	constraints = (
		sum(one_hot_constraints_layer) +
		sum(one_hot_constraints_num) +
		sum(CXgate_constraints)
	)

	cost = sum_poly(
		M-1,
		lambda m: sum_poly(
			N,
			lambda i: sum_poly(
				N,
				lambda j: sum_poly(
					N,
					lambda v: q[m][i][v]*q[m+1][j][v]
				) * ((N-1)*(i+j) - 2*i*j) / N
			)
		)
	)

	##########   solve   ##########

	solver = Solver(client)
	model = BinaryQuadraticModel(constraints * constraintWeight + cost)

	t1 = time.time()

	result = solver.solve(model)
	if len(result) == 0:
	    raise RuntimeError("Any one of constraints is not satisfied.")

	t2 = time.time()

	values = result[0].values
	q_values = decode_solution(q, values, 1)

	# print(q_values_main)

	##########   decode the result into string   ##########

	ans = []
	temp = []
	for m in range(M):
		s = ""
		t = [-1] * N
		for n in range(N):
			c = "."
			for v in range(N):
				if(q_values[m][n][v] > 0.5):
					c = chr(ord("a") + v)
					t[n] = v
			s += c
		ans.append(s)
		temp.append(t)

	cost = 0
	for m in range(M-1):
		cost += calcCost(temp[m], temp[m+1])

	t3 = time.time()

	# print("preparation time : " + str(t1-t0))
	# print("solving time : " + str(t2-t1))
	# print("decoding time : " + str(t3-t2))

	dt0, dt1, dt2 = t1-t0, t2-t1, t3-t2

	return dt0, dt1, dt2,cost, ans


# N = 3	# number of qubits
M = 20	# number of CX-gate layers

# mode = "STANDARD"
# mode = "Q_TIME"
mode = "COST"

filename = {
	"STANDARD" : "result_standard.txt",
	"Q_TIME" : "result_q_time.txt",
	"COST" : "result_cost.txt",
}
N_range = {
	"STANDARD" : range(3,27),
	"Q_TIME" : range(5,21,5),
	"COST" : range(3,7),
}

with open(filename[mode], mode = "w") as f:
	for N in N_range[mode]:
		print("N = " + str(N))
		f.write("N = " + str(N) + "\n")

		if(mode in ["Q_TIME", "COST"]):
			K = 10
		else:
			K = 1

		T0, T1, T2, Tsum = np.zeros(K), np.zeros(K), np.zeros(K), np.zeros(K)
		CC, CQ = np.zeros(K), np.zeros(K)

		for k in range(K):
			query = gen_CXgate_layers(N,M)
		
			if(mode == "STANDARD"):
				print(query)
				f.write(str(query) + "\n")

				print("")
				f.write("\n")
			
			else:
				print("k = " + str(k))


			if(N <= 7 and mode in ["STANDARD", "COST"]):
				dt, cost, ans = classic_solver(N,M,query)
				CC[k] = cost

				if(mode == "STANDARD"):
					print(ans)
					f.write(str(ans) + "\n")
					print("cost(classic) = " + str(cost))
					f.write("cost(classic) = " + str(cost) + "\n")
					print("time(classic) = " + str(dt))
					f.write("time(classic) = " + str(dt) + "\n")

					print("")
					f.write("\n")

			if(mode in ["STANDARD", "Q_TIME", "COST"]):
				dt0, dt1, dt2, cost, ans = quantum_solver(N,M,query)
				CQ[k] = cost

				T0[k] = dt0
				T1[k] = dt1
				T2[k] = dt2
				Tsum[k] = dt0+dt1+dt2

				if(mode == "STANDARD"):
					print(ans)
					f.write(str(ans) + "\n")
					print("cost(quantum) = " + str(cost))
					f.write("cost(quantum) = " + str(cost) + "\n")
					print("time(quantum) = " + str(Tsum[k]))
					f.write("time(quantum) = " + str(Tsum[k]) + "\n")

					print("")
					f.write("\n")

		if(mode == "Q_TIME"):
			print("preparation time : " + str(np.mean(T0)) + " +- " + str(np.std(T0, ddof=1)))
			print("solving time : " + str(np.mean(T1)) + " +- " + str(np.std(T1, ddof=1)))
			print("decoding time : " + str(np.mean(T2)) + " +- " + str(np.std(T2, ddof=1)))
			print("time(quantum) : " + str(np.mean(Tsum)) + " +- " + str(np.std(Tsum, ddof=1)))

			f.write("preparation time : " + str(np.mean(T0)) + " +- " + str(np.std(T0, ddof=1)) + "\n")
			f.write("solving time : " + str(np.mean(T1)) + " +- " + str(np.std(T1, ddof=1)) + "\n")
			f.write("decoding time : " + str(np.mean(T2)) + " +- " + str(np.std(T2, ddof=1)) + "\n")
			f.write("time(quantum) : " + str(np.mean(Tsum)) + " +- " + str(np.std(Tsum, ddof=1)) + "\n")

			print("")
			f.write("\n")

		if(mode == "COST"):
			print("cost(classic) : " + str(np.mean(CC)))
			print("cost(quantum) : " + str(np.mean(CQ)))

			f.write("cost(classic) : " + str(np.mean(CC)) + "\n")
			f.write("cost(quantum) : " + str(np.mean(CQ)) + "\n")

			print("")
			f.write("\n")







