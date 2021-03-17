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
client.parameters.timeout = 1000

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
			v0, v1 = query[m][g0], query[m][g0 + 1]

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

	result = solver.solve(model)
	if len(result) == 0:
		raise RuntimeError("Any one of constraints is not satisfied.")

	values = result[0].values
	q_values = decode_solution(q, values, 1)

	# print(q_values_main)

	##########   decode the result into string   ##########

	ans = [[-1 for n in range(N)] for m in range(M)]
	for m in range(M):
		for n in range(N):
			for v in range(N):
				if(q_values[m][n][v] > 0.5):
					ans[m][n] = v

	cost = 0
	for m in range(M-1):
		cost += calcCost(ans[m], ans[m+1])

	return cost, ans



##########   Main Application   ##########

import tkinter as tk
from tkinter import *
from tkinter import ttk

# Window Layout
WINDOW_WIDTH = 800
WINDOW_HEIGHT = 600

CANVAS_WIDTH = WINDOW_WIDTH * 3
CANVAS_HEIGHT = WINDOW_HEIGHT * 0.6
CANVAS_Y = WINDOW_HEIGHT * 0.25

label_head = WINDOW_WIDTH * 0.15
comboBox_head = WINDOW_WIDTH * 0.35
button_head = WINDOW_WIDTH * 0.7

row_head_1 = WINDOW_HEIGHT * 0.05
row_head_2 = WINDOW_HEIGHT * 0.10
row_head_3 = WINDOW_HEIGHT * 0.15

# Variables
canvas = None
qasmFileName_input = None
timeout_input = None
constraintWeight_input = None

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
def convertToColorCode(r, g, b):
	r = max(min(r, 255), 0)
	g = max(min(g, 255), 0)
	b = max(min(b, 255), 0)

	return "#" + format(int(r), '02x') +  format(int(g), '02x') + format(int(b), '02x')

def getColor(theta_deg):
	while(theta_deg < 0): theta_deg += 360
	while(theta_deg > 360): theta_deg -= 360

	r,g,b = 0,0,0

	if(theta_deg < 60):
		r = 255
		g = 255 * theta_deg / 60
		b = 0
	elif(theta_deg < 120):
		r = 255 * (120 - theta_deg) / 60
		g = 255
		b = 0
	elif(theta_deg < 180):
		r = 0
		g = 255
		b = 255 * (theta_deg - 120) / 60
	elif(theta_deg < 240):
		r = 0
		g = 255 * (240 - theta_deg) / 60
		b = 255
	elif(theta_deg < 300):
		r = 255 * (theta_deg - 240) / 60
		g = 0
		b = 255
	else:
		r = 255
		g = 0
		b = 255 * (360 - theta_deg) / 60

	return convertToColorCode(r,g,b)

# Canvas Layout
CANVAS_X_START = 50
CANVAS_X_OFFSET = 50

DRAW_LAYER_PARTITION = True

def Draw(_gates):
	initial_symbols = _gates[0]
	_gates = _gates[1:]


	##########   gates -> draw_layers   ##########

	gate_block = [[]]
	for gate in _gates:
		if(gate[0] == "partition"):
			gate_block.append([])
		else:
			gate_block[-1].append(gate)

	# print(gate_block)

	draw_layers = []
	symbols = [s for s in initial_symbols]

	for m in range(len(gate_block)):
		gates = gate_block[m]
		Ng = len(gates)
		used = [False] * Ng

		for _i in range(Ng):
			# 0 -> empty
			# 1 -> reserved
			# 2 -> used
			stat = [0] * N
			layer = []

			for i in range(Ng):
				if(used[i] == True): continue
				g = gates[i]

				if(g[0] == "cx" or g[0] == "swap"):
					flag = True
					j0, j1 = symbols.index(g[1]), symbols.index(g[2])
					
					if(stat[j0] == 1 or stat[j1] == 1):
						flag = False
					for j in range(min(j0,j1), max(j0,j1)+1):
						if(stat[j] == 2):
							flag = False
							break

					if(flag):
						used[i] = True
						layer.append(g)
						for j in range(min(j0,j1), max(j0,j1)+1):
							stat[j] = 2
					else:
						stat[j0] = max(stat[j0], 1)
						stat[j1] = max(stat[j1], 1)

				elif(g[0] == "u3"):
					j = symbols.index(g[1])
					
					if(stat[j] == 0):
						used[i] = True
						layer.append(g)
						stat[j] = 2

			if(layer == []): break

			for g in layer:
				if(g[0] == "swap"):
					i, j = symbols.index(g[1]), symbols.index(g[2])
					symbols[i],symbols[j] = symbols[j],symbols[i]

			draw_layers.append(layer)

		if(DRAW_LAYER_PARTITION == True):
			draw_layers.append([["partition"]])

	if(DRAW_LAYER_PARTITION == True):
		draw_layers = draw_layers[:-1]
	
	# print(draw_layers)


	##########   draw   ##########

	GATE_SIZE = 10

	Y = [(CANVAS_HEIGHT-20) * (n+0.5) / N for n in range(N)]
	colors = [getColor(360 * n/N) for n in range(N)]
	CX_COLOR = getColor(190)
	U3_COLOR = getColor(310)

	canvas.delete("all")

	for n in range(N):
		canvas.create_text(
			CANVAS_X_START - 10, Y[n], text = "q"+str(initial_symbols[n]),
			font = ("", 20), anchor = "e", fill = colors[initial_symbols[n]]
		)
		canvas.create_line(
			CANVAS_X_START, Y[n], CANVAS_X_START + CANVAS_X_OFFSET, Y[n],
			fill = colors[initial_symbols[n]], width = 3
		)

	Nl = len(draw_layers)
	prev_symbols = initial_symbols

	for n in range(Nl):
		xd = (CANVAS_WIDTH - CANVAS_X_START - CANVAS_X_OFFSET * 2) / (Nl-1)
		x = CANVAS_X_START + CANVAS_X_OFFSET + n * xd
		xl, xr = x - xd/2, x + xd/2

		current_symbols = [s for s in prev_symbols]
		for g in draw_layers[n]:
			if(g[0] == "swap"):
				i,j = current_symbols.index(g[1]),current_symbols.index(g[2])
				current_symbols[i],current_symbols[j] = current_symbols[j],current_symbols[i]

		for i in range(N):
			canvas.create_line(
				xl,Y[i],x,Y[i],
				fill = colors[prev_symbols[i]], width = 3
			)

		for i in range(N):
			canvas.create_line(
				x,Y[i],xr,Y[i],
				fill = colors[current_symbols[i]], width = 3
			)
		
		for g in draw_layers[n]:
			if(g[0] == "u3"):
				y = Y[prev_symbols.index(g[1])]
				canvas.create_rectangle(
					x-GATE_SIZE, y-GATE_SIZE, x+GATE_SIZE, y+GATE_SIZE,
					fill = U3_COLOR
				)

			elif(g[0] == "cx"):
				y1,y2 = Y[prev_symbols.index(g[1])],Y[prev_symbols.index(g[2])]
				canvas.create_line(
					x,y1,x,y2,
					fill = CX_COLOR, width = 2
				)
				canvas.create_oval(
					x-GATE_SIZE/2, y1-GATE_SIZE/2, x+GATE_SIZE/2, y1+GATE_SIZE/2,
					fill = CX_COLOR
				)
				canvas.create_oval(
					x-GATE_SIZE, y2-GATE_SIZE, x+GATE_SIZE, y2+GATE_SIZE,
					fill = CX_COLOR
				)
				canvas.create_line(
					x-GATE_SIZE*0.7,y2,x+GATE_SIZE*0.7,y2,
					fill = "white", width = 3
				)
				canvas.create_line(
					x,y2-GATE_SIZE*0.7,x,y2+GATE_SIZE*0.7,
					fill = "white", width = 3
				)

			elif(g[0] == "swap"):
				y1,y2 = Y[prev_symbols.index(g[1])],Y[prev_symbols.index(g[2])]
				canvas.create_line(
					x,y1,x,y2,
					fill = CX_COLOR, width = 2
				)
				canvas.create_line(
					x-GATE_SIZE*0.7, y2-GATE_SIZE*0.7, x+GATE_SIZE*0.7, y2+GATE_SIZE*0.7,
					fill = CX_COLOR, width = 2
				)
				canvas.create_line(
					x+GATE_SIZE*0.7, y2-GATE_SIZE*0.7, x-GATE_SIZE*0.7, y2+GATE_SIZE*0.7,
					fill = CX_COLOR, width = 2
				)
				canvas.create_line(
					x-GATE_SIZE*0.7, y1-GATE_SIZE*0.7, x+GATE_SIZE*0.7, y1+GATE_SIZE*0.7,
					fill = CX_COLOR, width = 2
				)
				canvas.create_line(
					x+GATE_SIZE*0.7, y1-GATE_SIZE*0.7, x-GATE_SIZE*0.7, y1+GATE_SIZE*0.7,
					fill = CX_COLOR, width = 2
				)

			elif(g[0] == "partition"):
				canvas.create_line(
					x, 0, x, CANVAS_HEIGHT, dash = (10, 5)
				)

		prev_symbols = current_symbols


def Load():
	qasmFileName = qasmFileName_input.get()

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
			Q_name, Q_symbol = qv[0]+"["+str(i)+"]", n+i
			Q_dict[Q_name] = Q_symbol
			Q_dict_inv[Q_symbol] = Q_name

	global N
	N = len(Q_dict)

	# print(Q_dict)
	# print(Q_dict_inv)

	##########   detect gates   ##########
	global gates_input
	gates_input = []

	gates_input.append([])
	for n in range(N):
		gates_input[0].append(n)

	for s in S:
		s = s.replace(","," ")
		s = s.replace(";"," ")
		for Q_name, Q_symbol in Q_dict.items():
			s = s.replace(Q_name,str(Q_symbol))

		A = s.split()
		if(A[0] == "cx"): gates_input.append(["cx", int(A[1]), int(A[2])])
		elif(A[0][:2] == "u3"): gates_input.append(["u3", int(A[-1])])

	#print(gates_input)
	#Draw(gates_input)


	##########   detect layers   ##########

	global gates_layer
	gates_layer = [gates_input[0]]

	Ng = len(gates_input) - 1
	used = [False] * Ng

	for _i in range(Ng):
		flag = False
		stat = [False] * N

		for i in range(Ng):
			if(used[i] == True): continue
			g = gates_input[i+1]

			if(g[0] == "cx" or g[0] == "swap"):
				j0, j1 = g[1], g[2]
				
				if(stat[j0] == False and stat[j1] == False):
					used[i] = True
					flag = True
					gates_layer.append(g)

				stat[j0] = True
				stat[j1] = True

			else:
				j = g[1]
				if(stat[j] == False):
					used[i] = True
					flag = True
					gates_layer.append(g)

		if(flag == False): break

		gates_layer.append(["partition"])

	gates_layer = gates_layer[:-1]

	for i in range(len(gates_layer))[::-1]:
		layer = gates_layer[i]
		if(layer[0] == "cx"):
			break
		elif(layer[0] == "partition"):
			gates_layer.pop(i)
			break

	# print(gates_layer)
	Draw(gates_layer)

def bubble_sort(symbols_0, symbols_1):
	s0 = [s for s in symbols_0]
	s1 = [s for s in symbols_1]
	ans = []
	n = len(s0)

	for j in range(n)[::-1]:
		j0 = s0.index(s1[j])

		for i in range(j0, j):
			ans.append(["swap", s0[i], s0[i+1]])
			s0[i],s0[i+1] = s0[i+1],s0[i]

	return ans

def Execute():
	if(gates_layer == []):
		return

	global client, constraintWeight
	client.parameters.timeout = int(timeout_input.get().replace("ms",""))
	constraintWeight = int(constraintWeight_input.get())

	global M,query
	query = [[]]

	for layer in gates_layer[1:]:
		if(layer[0] == "partition"):
			query.append([])
		elif(layer[0] == "cx"):
			query[-1].append(layer[1])
			query[-1].append(layer[2])

	M = len(query)

	# print(N)
	# print(M)
	# print(query)

	cost, ans = quantum_solver(N,M,query)

	print("cost = " + str(cost))
	# print(ans)

	global gates_swap
	gates_swap = []

	gates_swap.append(ans[0])
	m = 0
	for layer in gates_layer[1:]:
		gates_swap.append(layer)

		if(layer[0] == "partition"):
			m += 1
			temp = bubble_sort(ans[m-1], ans[m])
			for t in temp:
				gates_swap.append(t)

			gates_swap.append(["partition"])

	# print(gates_swap)
	Draw(gates_swap)

# Viewer
root = tk.Tk()
root.title("QC Viewer for LNNA")
root.geometry(str(WINDOW_WIDTH) + "x" + str(WINDOW_HEIGHT))
root.resizable(False, False)

canvas = tk.Canvas(root, bg = "white")
canvas.place(x = 0, y = CANVAS_Y, width = WINDOW_WIDTH, height = CANVAS_HEIGHT)
bar_x = tk.Scrollbar(canvas, orient = tk.HORIZONTAL)
bar_x.pack(side = tk.BOTTOM, fill = tk.X)
bar_x.config(command = canvas.xview)
canvas.config(xscrollcommand = bar_x.set)
canvas.config(scrollregion = (0, 0, CANVAS_WIDTH, WINDOW_HEIGHT))

# QASM File Input
qasmFileName_input = StringVar()
qasmFileName_height = row_head_1
qasmFileName_label = tk.Label(text = 'サンプル名')
qasmFileName_label.place(x = label_head, y = qasmFileName_height)

qasmFileName_comboBox = ttk.Combobox(root, textvariable = qasmFileName_input)
qasmFileName_comboBox['values'] = ('qasm/ex1.txt','qasm/ex2.txt','qasm/ex3.txt')
qasmFileName_comboBox.set("qasm/ex1.txt")
qasmFileName_comboBox.place(x = comboBox_head, y = qasmFileName_height)

# Timeout Input
timeout_input = StringVar()
timeout_height = row_head_2
timeout_label = tk.Label(text = 'タイムアウト')
timeout_label.place(x = label_head, y = timeout_height)

timeout_comboBox = ttk.Combobox(root, textvariable = timeout_input)
timeout_comboBox['values'] = ('1000ms','2000ms','3000ms','5000ms','10000ms','20000ms','30000ms')
timeout_comboBox.set(str(client.parameters.timeout) + "ms")
timeout_comboBox.place(x = comboBox_head, y = timeout_height)

# Constraint Weight Input
constraintWeight_input = StringVar()
constraintWeight_height = row_head_3
constraintWeight_label = tk.Label(text = '制約の重みパラメータ')
constraintWeight_label.place(x = label_head, y = constraintWeight_height)

constraintWeight_comboBox = ttk.Combobox(root, textvariable = constraintWeight_input)
constraintWeight_comboBox['values'] = ('1','10','100','1000','10000')
constraintWeight_comboBox.set(str(constraintWeight))
constraintWeight_comboBox.place(x = comboBox_head, y = constraintWeight_height)

# Load Button
loadButton = tk.Button(root, text = '読み込む', command = Load)
loadButton.place(x = button_head, y = row_head_1 - 5, width = 100)

# Execute Button
executeButton = tk.Button(root, text = '実行', command = Execute)
executeButton.place(x = button_head, y = row_head_3 - 5, width = 100)

# Main
root.mainloop()



