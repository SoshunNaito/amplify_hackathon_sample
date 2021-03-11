import itertools
import matplotlib.pyplot as plt

n = 8

def calcCost(A):
	ans = 0
	for i in range(n):
		for j in range(i+1,n):
			if(A[i]>A[j]):
				ans += 1
	return ans

def func(i,j):
	return (i+j)/2 - i*j/(n-1)

def test(A):
	ans = 0
	for i in range(n):
		ans += func(i, A[i])
	return ans


X,Y = [],[]

for p in itertools.permutations(range(n), n):
	X.append(calcCost(p))
	Y.append(test(p))

plt.scatter(X,Y)
plt.xlim(0,n*(n-1)//2)
plt.ylim(0,n*(n-1)//2)
plt.show()