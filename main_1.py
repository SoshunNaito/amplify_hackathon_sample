import numpy as np
import random
import itertools
import queue

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

client = FixstarsClient()
# client.token = "DELETED TOKEN"
client.token = "Xiccn8dKHhDoboWnaixrUEDRjvMl2vzo"
client.parameters.timeout = 10 * 1000

N = 4	# number of qubits
M = 2	# number of CX-gate layers

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

def classic_solver(query): # solve with a classic computer
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
	dp = [[INF for f in range(F)] for m in range(M)]
	back = [[-1 for f in range(F)] for m in range(M)]

	for m in range(M):
		if(m == 0):
			for f in range(F):
				dp[m][f] = 0

				valid = True
				for g0 in range(0,len(query[m]),2):
					v0, v1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f].find(v0)-x[f].find(v1)) != 1):
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
					v0, v1 = query[m][g0], query[m][g0 + 1]
					if(abs(x[f2].find(v0)-x[f2].find(v1)) != 1):
						valid = False
						break

				if(valid == False):
					dp[m][f2] = INF
					continue

				##########   01-BFS   ##########
				dist = [-1 for f in range(F)]
				dist[f2] = 0

				que = queue.Queue()
				que.put(f2)

				while not que.empty():
					f1 = que.get()
					s = x[f1]
					for i in range(N-1):
						v0, v1 = s[i], s[i+1]
						s = s.translate(str.maketrans({v0:v1, v1:v0}))

						f = x_inv[s]
						if(dist[f] == -1):
							dist[f] = dist[f1] + 1
							que.put(f)
							if(dp[m-1][f] + dist[f] < dp[m][f2]):
								back[m][f2] = f
								dp[m][f2] = dp[m-1][f] + dist[f]

						s = s.translate(str.maketrans({v0:v1, v1:v0}))

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
	return cost, ans

def quantum_solver(query): # solve with Amplify
	q_main = gen_symbols(BinaryPoly, M, N, N)	# represent the solution
	q_sub = gen_symbols(BinaryPoly, M-1, N, N)	# ancilla qubits

	##########   constraints   ##########

	# each layer doesn't have 2+ same values
	one_hot_constraints_layer_main = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q_main[m][n][v] for n in range(N)), 1)
		for m in range(M)
		for v in range(N)
	]

	# each qubit doesn't have 2+ values
	one_hot_constraints_num_main = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q_main[m][n][v] for v in range(N)), 1)
		for m in range(M)
		for n in range(N)
	]

	# Every CX-gate must be applied for 2 adjacent qubits
	CXgate_constraints_main = []
	for m in range(M):
		for g0 in range(0,len(query[m]),2):
			v0, v1 = ord(query[m][g0])-ord("a"), ord(query[m][g0 + 1])-ord("a")

			# v0 and v1 must be adjacent each other
			for i in range(N):
				for j in range(i+2, N):
					CXgate_constraints_main.append(
						penalty(q_main[m][i][v0] * q_main[m][j][v1])
					)
					CXgate_constraints_main.append(
						penalty(q_main[m][i][v1] * q_main[m][j][v0])
					)

	# each column doesn't have 2+ values
	one_hot_constraints_column_sub = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q_sub[m][n][v] for n in range(N)), 1)
		for m in range(M-1)
		for v in range(N)
	]

	# each row doesn't have 2+ values
	one_hot_constraints_row_sub = [
		# m -> layer
		# n -> qubit
		# v -> value of qubit
		equal_to(sum(q_sub[m][n][v] for v in range(N)), 1)
		for m in range(M-1)
		for n in range(N)
	]

	# q_sub qubits represent "C" in calcCost(A,B)
	matrix_constraints_sub = [
		equal_to(sum([q_main[m][i][c] * q_main[m+1][j][c] for c in range(N)]) - q_sub[m][i][j], 0)
		for i in range(N)
		for j in range(N)
		for m in range(M-1)
	]

	constraints = (
		sum(one_hot_constraints_layer_main) +
		sum(one_hot_constraints_num_main) +
		sum(CXgate_constraints_main) +
		sum(one_hot_constraints_column_sub) +
		sum(one_hot_constraints_row_sub) +
		sum(matrix_constraints_sub)
	)

	##########   solve   ##########

	solver = Solver(client)

	model = BinaryQuadraticModel(constraints)
	result = solver.solve(model)
	if len(result) == 0:
	    raise RuntimeError("Any one of constraints is not satisfied.")

	values = result[0].values
	q_values_main = decode_solution(q_main, values, 1)
	q_values_sub = decode_solution(q_sub, values, 1)

	# print(q_values_main)

	##########   decode the result into string   ##########

	ans = []
	for m in range(M):
		s = ""
		for n in range(N):
			c = "."
			for v in range(N):
				if(q_values_main[m][n][v] > 0.5):
					c = chr(ord("a") + v)
			s += c
		ans.append(s)

	# print(ans)

	ans_sub = []
	cost = 0
	for m in range(M-1):
		s = ""
		for n in range(N):
			c = "."
			for v in range(N):
				if(q_values_sub[m][n][v] > 0.5):
					c = str(v)
			s += c
		for i in range(N):
			for j in range(i+1,N):
				if(ord(s[i]) > ord(s[j])):
					cost += 1
		ans_sub.append(s)

	# print(ans)
	# print(ans_sub)
	
	return cost, ans

RUN_CLASSIC_SOLVER = True
RUN_QUANTUM_SOLVER = True

query = gen_CXgate_layers()
print(query)

if(RUN_CLASSIC_SOLVER):
	cost, ans = classic_solver(query)
	print("cost(classic) = " + str(cost))
	print(ans)

if(RUN_QUANTUM_SOLVER):
	timeout = 60 / 32
	timeout_max = 30 + 0.0001

	while(True):
		client.parameters.timeout = int(timeout * 1000 + 0.5)

		try:
			print("timeout = " + str(timeout) + " secs")
			cost, ans = quantum_solver(query)
			print("cost(quantum) = " + str(cost))
			print(ans)
			break

		except RuntimeError as e:
			timeout *= 2
			if(timeout > timeout_max):
				break







