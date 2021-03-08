import numpy as np
from math import sqrt

from amplify import Solver, decode_solution, gen_symbols, BinaryPoly, sum_poly, BinaryQuadraticModel
from amplify.client import FixstarsClient
from amplify.constraint import equal_to, penalty

client = FixstarsClient()
client.token = "Xiccn8dKHhDoboWnaixrUEDRjvMl2vzo"
client.parameters.timeout = 1000


sudoku_initial = np.array(
    [
        [2, 0, 5, 1, 3, 0, 0, 0, 4],
        [0, 0, 0, 0, 4, 8, 0, 0, 0],
        [0, 0, 0, 0, 0, 7, 0, 2, 0],
        [0, 3, 8, 5, 0, 0, 0, 9, 2],
        [0, 0, 0, 0, 9, 0, 7, 0, 0],
        [0, 0, 0, 0, 0, 0, 4, 5, 0],
        [8, 6, 0, 9, 7, 0, 0, 0, 0],
        [9, 5, 0, 0, 0, 0, 0, 3, 1],
        [0, 0, 4, 0, 0, 0, 0, 0, 0],
    ]
)

def print_sudoku(sudoku):
    
    N = len(sudoku)
    n = int(sqrt(N))
    width = len(str(N))
    for i in range(N):
        line = ""
        if i % n == 0 and i != 0:
            print("-" * ((width + 1) * n * n + 2 * (n - 1) - 1))
        for j in range(len(sudoku[i])):
            if j % n == 0 and j != 0:
                line += "| "
            line += (str(sudoku[i][j]).rjust(width) if sudoku[i][j] != 0 else " ") + " "
        print(line)

def sudoku_solve(*args):

    N = len(sudoku_initial)
    n = int(sqrt(N))
    q = gen_symbols(BinaryPoly, N, N, N)

    for i, j in zip(*np.where(sudoku_initial != 0)):
        k = sudoku_initial[i, j] - 1

        for m in range(N):
            q[i][m][k] = BinaryPoly(0)
            q[m][j][k] = BinaryPoly(0)
            q[i][j][m] = BinaryPoly(0)
            q[n * (i // n) + m // n][n * (j // n) + m % n][k] = BinaryPoly(0)  # 制約(c)

        q[i][j][k] = BinaryPoly(1)

    row_constraints = [
        equal_to(sum(q[i][j][k] for j in range(N)), 1) for i in range(N) for k in range(N)
    ]

    col_constraints = [
        equal_to(sum(q[i][j][k] for i in range(N)), 1) for j in range(N) for k in range(N)
    ]

    num_constraints = [
        equal_to(sum(q[i][j][k] for k in range(N)), 1) for i in range(N) for j in range(N)
    ]

    block_constraints = [
        equal_to(sum([q[i + m // n][j + m % n][k] for m in range(N)]), 1)
        for i in range(0, N, n)
        for j in range(0, N, n)
        for k in range(N)
    ]

    constraints = (
        sum(row_constraints)
        + sum(col_constraints)
        + sum(num_constraints)
        + sum(block_constraints)
    )

    solver = Solver(client)
    solver.client.parameters.timeout = 1000
    model = BinaryQuadraticModel(constraints)
    result = solver.solve(model)
    i = 0
    while len(result) == 0:
        if i > 5:
            raise RuntimeError()
        result = solver.solve(model)
        i += 1

    values = result[0].values
    q_values = decode_solution(q, values)
    answer = np.array([np.where(np.array(q_values[i]) != 0)[1] + 1 for i in range(N)])
    print_sudoku(answer)

print("initial")
print_sudoku(sudoku_initial)
print("solved")
sudoku_solve(sudoku_initial)