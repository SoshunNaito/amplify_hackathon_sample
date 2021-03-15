import numpy as np
import random
import itertools
import queue

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

import time

##########   Amplify Subroutine   ##########

client = FixstarsClient()
# client.token = "DELETED TOKEN"
client.token = "Xiccn8dKHhDoboWnaixrUEDRjvMl2vzo"
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



##########   Main Application   ##########

import tkinter as tk
from tkinter import *
from tkinter import ttk

WINDOW_WIDTH = 800
WINDOW_HEIGHT = 600

root = tk.Tk()

# Variables
canvas = None
qasmFileName_input = StringVar()
timeout_input = StringVar()
constraintWeight_input = StringVar()

N = 10
M = 10
query = []

Q_variables = []
Q_dict = {}
Q_dict_inv = {}

gates_input = []
gates_layer = []
gates_swap = []

# Functions
def getColor(theta_deg):
	while(theta_deg < 0): theta_deg += 360
	while(theta_deg > 360): theta_deg -= 360

	ans = "#"
	if(theta_deg < 60):
		ans += "ff"
		ans += format((int)(255 * theta_deg / 60), '02x')
		ans += "00"
	elif(theta_deg < 120):
		ans += format((int)(255 * (120 - theta_deg) / 60), '02x')
		ans += "ff"
		ans += "00"
	elif(theta_deg < 180):
		ans += "00"
		ans += "ff"
		ans += format((int)(255 * (theta_deg - 120) / 60), '02x')
	elif(theta_deg < 240):
		ans += "00"
		ans += format((int)(255 * (240 - theta_deg) / 60), '02x')
		ans += "ff"
	elif(theta_deg < 300):
		ans += format((int)(255 * (theta_deg - 240) / 60), '02x')
		ans += "00"
		ans += "ff"
	else:
		ans += "ff"
		ans += "00"
		ans += format((int)(255 * (360 - theta_deg) / 60), '02x')

	return ans

def Draw(gates):
	CANVAS_WIDTH = canvas.winfo_width()
	CANVAS_HEIGHT = canvas.winfo_height() - 20

	symbols = gates[0]
	gates = gates[1:]

	X = 30
	Y = [CANVAS_HEIGHT * (n+0.5) / N for n in range(N)]
	colors = [getColor(360 * n/N) for n in range(N)]

	for n in range(N):
		canvas.create_text(X - 10, Y[n], text = symbols[n], font = ("", 20), anchor = "e", fill = colors[n])
		canvas.create_line(X, Y[n], X + 100, Y[n], fill = colors[n], width = 3)

def Load():
	qasmFileName = qasmFileName_input.get() + ".txt"

	S = []
	with open(qasmFileName, mode = "r") as f:
		S = f.readlines()

	##########   detect variables   ##########
	global Q_variables, Q_dict, Q_dict_inv
	Q_variables = []
	Q_dict = {}
	Q_dict_inv = {}

	for s in S:
		A = s.split()
		if(A[0]=="qreg"):
			i,j = A[1].find("["), A[1].find("]")
			Q_variables.append([A[1][:i], int(A[1][i+1:j])])

	# print(Q_variables)
	
	for qv in Q_variables:
		n = len(Q_dict)
		for i in range(qv[1]):
			Q_name, Q_symbol = qv[0]+"["+str(i)+"]", chr(ord("a")+n+i)
			Q_dict[Q_name] = Q_symbol
			Q_dict_inv[Q_symbol] = Q_name

	global N
	N = len(Q_dict)

	# print(Q_dict)
	# print(Q_dict_inv)

	##########   detect CX-gates   ##########
	global gates_input
	gates_input = []

	gates_input.append([])
	for n in range(N):
		gates_input[0].append(chr(ord("a")+n))

	for s in S:
		s = s.replace(","," ")
		s = s.replace(";"," ")
		for Q_name, Q_symbol in Q_dict.items():
			s = s.replace(Q_name,Q_symbol)

		A = s.split()
		if(A[0] == "cx"): gates_input.append(["cx", A[1], A[2]])
		elif(A[0][:2] == "u3"): gates_input.append(["u3", A[-1]])

	#print(gates_input)
	Draw(gates_input)

def Execute():
	print(gates_input)

# Window Layout
root.title("QC Viewer for LNNA")
root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT))
root.resizable(False, False)

label_head = WINDOW_WIDTH * 0.15
comboBox_head = WINDOW_WIDTH * 0.35
button_head = WINDOW_WIDTH * 0.7

row_head_1 = WINDOW_HEIGHT * 0.05
row_head_2 = WINDOW_HEIGHT * 0.10
row_head_3 = WINDOW_HEIGHT * 0.15

# Viewer
canvas = tk.Canvas(root, bg = "#bbb")
canvas.place(x = 0, y = WINDOW_HEIGHT * 0.25, width = WINDOW_WIDTH, height = WINDOW_HEIGHT * 0.5)
bar_x = tk.Scrollbar(canvas, orient = tk.HORIZONTAL)
bar_x.pack(side = tk.BOTTOM, fill = tk.X)
bar_x.config(command = canvas.xview)
canvas.config(xscrollcommand = bar_x.set)
canvas.config(scrollregion = (0, 0, WINDOW_WIDTH * 2, WINDOW_HEIGHT))

# QASM File Input
qasmFileName_height = row_head_1
qasmFileName_label = tk.Label(text = 'サンプル名')
qasmFileName_label.place(x = label_head, y = qasmFileName_height)

qasmFileName_comboBox = ttk.Combobox(root, textvariable = qasmFileName_input)
qasmFileName_comboBox['values'] = ('qasm/ex1','qasm/ex2','qasm/ex3')
qasmFileName_comboBox.set("qasm/ex1")
qasmFileName_comboBox.place(x = comboBox_head, y = qasmFileName_height)

# Timeout Input
timeout_height = row_head_2
timeout_label = tk.Label(text = 'タイムアウト')
timeout_label.place(x = label_head, y = timeout_height)

timeout_comboBox = ttk.Combobox(root, textvariable = timeout_input)
timeout_comboBox['values'] = ('1000ms','2000ms','3000ms','5000ms','10000ms','20000ms','30000ms')
timeout_comboBox.set("1000ms")
timeout_comboBox.place(x = comboBox_head, y = timeout_height)

# Constraint Weight Input
constraintWeight_height = row_head_3
constraintWeight_label = tk.Label(text = '制約の重みパラメータ')
constraintWeight_label.place(x = label_head, y = constraintWeight_height)

constraintWeight_comboBox = ttk.Combobox(root, textvariable = constraintWeight_input)
constraintWeight_comboBox['values'] = ('x1','x10','x100','x1000','x10000')
constraintWeight_comboBox.set("x100")
constraintWeight_comboBox.place(x = comboBox_head, y = constraintWeight_height)

# Load Button
loadButton = tk.Button(root, text = '読み込む', command = Load)
loadButton.place(x = button_head, y = row_head_1 - 5, width = 100)

# Execute Button
executeButton = tk.Button(root, text = '実行', command = Execute)
executeButton.place(x = button_head, y = row_head_3 - 5, width = 100)

# Main
root.mainloop()



