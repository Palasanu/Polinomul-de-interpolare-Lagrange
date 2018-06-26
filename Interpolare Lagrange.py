from numpy import matrix,array,allclose,dot,matmul
from numpy.linalg import norm, solve, inv
from copy import deepcopy
import random

x = [0, 1, 2, 3, 4, 5]
f = [50, 47, -2, -121, -310, -545]

x_prim = 1.5

def calc_h(n,x0,xn):
	return (xn-x0)/n

def calc_t(x_prim,x,h):
	return (x_prim-x[0])/h

def func(x):
	return (x**4) - (12 * (x**3)) + (30 * (x**2)) + 12

def pasulk(f):
	n = len(f)
	y = list()
	for i in range(n-1):
		y.append(f[i+1]-f[i])
	return y

def SchemaAitken(f):
	n = len(f)
	y = list()
	yk = f.copy()
	y.append(f[0])
	print("Aitken, pasul 0 : ",yk )
	for i in range(n-1):
		yk = pasulk(yk)
		print("Aitken, pasul",i+1,": ",yk )
		y.append(yk[0])
	print("y0 in Schema lui Aitken: ",y)
	return y

def calc_sk(t,n):
	sk = [None,t]
	for k in range(2,n):
		aux = sk[k-1]*((t-k+1)/k)
		sk.append(aux)
	return sk

def afisare_matrice(m, pading=0, spaces = 0):		#matrice e o lista de liste, [int] padding spatiul de la stanga, [int]spaces distanta dintre coloane
	elements = list()
	matrice = deepcopy(m)
	for linie in matrice:
		for el in linie:
			el = round(el,2)
			elements.append(len(str(el)))
	l = max(elements)
	for linie in matrice:
		print(' '*pading, end='')
		for el in linie:
			str_el = str(round(el,2))
			pad = l - len(str_el)
			print(' '*pad + str_el + ' '*spaces, end=' ')
		print()

def make_Ab(A,b): #concateaza matricea A cu matricea b
	n = len(A)
	Ab = deepcopy(A)
	for index in range(0,n):
		Ab[index].append(b[index][0])
	return deepcopy(Ab)


def pivotare(M,l): #cauta pivotul partial
	n = len(M)
	A = deepcopy(M)
	i_max = l
	for i in range(l+1,n):
		if abs(A[i][l]) > abs(A[i_max][l]) :
			i_max = i
	if abs(A[l][l]) < abs(A[i_max][l]):
		A[l], A[i_max] = A[i_max].copy(), A[l].copy()
	return deepcopy(A)

def transformare(M,l): #transforma matricea in una echivalenta cu coloana de sub pivot in 0
	n = len(M)
	A = deepcopy(M)
	for i in range(l+1,n):
		aux = M[l].copy()	
		for j in range(l,n+1):
			aux[j] = aux[j] * (M[i][l] / M[l][l])
			A[i][j] = M[i][j] - aux[j]
	return A

def triunghiularizare(Ab): #transorma matricea in una triunchiulara sus
	n = len(Ab)
	for l in range(0,n): #pasul l
		Ab = pivotare(Ab,l)
		if Ab[l][l] == 0:
			return False
		Ab = transformare(Ab,l)
	return deepcopy(Ab)

def afisare_matrice(m, pading=0, spaces = 0):		#matrice e o lista de liste, [int] padding spatiul de la stanga, [int]spaces distanta dintre coloane
	elements = list()
	matrice = deepcopy(m)
	for linie in matrice:
		for el in linie:
			el = round(el,2)
			elements.append(len(str(el)))
	l = max(elements)
	for linie in matrice:
		print(' '*pading, end='')
		for el in linie:
			str_el = str(round(el,2))
			pad = l - len(str_el)
			print(' '*pad + str_el + ' '*spaces, end=' ')
		print()


def rezolva_sistem(Ab): #rezolva sistemul si retuneaza o lista cu solutii
	n = len(Ab)
	x = [0]*n
	x[n-1] = Ab[n-1][n] / Ab[n-1][n-1]
	for i in range(n-1,-1,-1):
		sum = 0
		for j in range(n-1,i,-1):
			sum += Ab[i][j] * x[j]
		x[i] = (Ab[i][n] - sum) / Ab[i][i]
	sol = list()
	for el in x:
		sol.append([el])
	return sol.copy()


def initializeaza(a):
	with open('initializare.txt') as f:
		line = f.readline()
		aux = line.split(' ')
	n = int(aux[0])
	x0 = int(aux[1])
	xn = int(aux[2])
	x = list()
	x.append(x0)
	f = list()
	h = calc_h(n,x0,xn)
	for i in range(1,xn-1):
		x.append(x[-1]+h)
	x.append(xn)
	for el in x:
		f.append(d(len(a)-1,a,el))

	return x,f

def interpolareLagrange(x_prim,x,f):
	print("  x :",x)
	print("f(x):",f)

	n = len(f)
	h = calc_h(n-1,x[0],x[-1])
	t = calc_t(x_prim,x,h)

	print("h: ", h)
	print("t: ", t)

	y = SchemaAitken(f)
	sk = calc_sk(t,n)

	print("sk: ",sk)

	rez = y[0]
	for i in range(1,n):
		rez += y[i]*sk[i]

	print("f(",x_prim,") =",rez)
	return rez

def metoda(x_prim,x,y,m):
	m=m+1
	B = list()
	for i in range(m):
		B.append([])
		for j in range(m):
			sum = 0
			for k in range(m):
				sum += x[k]**(i+j)	
			B[i].append(sum)

	f = list()
	for i in range(0,m):
		sum=0
		for k in range(0,m):
			sum = y[k]*x[k]**i
		f.append([sum])

	# M = array(B)
	# n = array(f)
	# a = solve(B,f)

	# print(B)
	# print()
	# print(f)

	Bf = make_Ab(B,f)
	afisare_matrice(Bf,2,2)
	Bf = triunghiularizare(Bf)
	a = rezolva_sistem(Bf)

	afiseazaPolinom(a)
	sol = d(len(a)-1,a,x_prim)
	print("Solutia f(x) =",sol)
	return a

def d(i,c,x):
	if i == 0:
		return c[0][0]
	else:
		return c[i][0] + (d(i-1,c,x)*x)

def Q(x,p,c):
	rez = d(p-1)
	for i in range(1,p):
		rez += d(i-1,c,x)*x**p-i

def afiseazaPolinom(x):
	print("f(x) = ",end='')
	n = len(x)
	for i in range(n,0,-1):
		if x[n-i][0] > 1:
			if i-1 == 1:
				print( " +" + str(round(x[n-i][0],3)) + "*x"  ,end = "")
			elif i-1 == 0:
				print( " +" + str(round(x[n-i][0],3)))
			elif i == n:
				print( str(round(x[n-i][0],3)) + "*x^" + str(i-1) ,end = "")
			else:
				print( " +" + str(round(x[n-i][0],3)) + "*x^" + str(i-1) ,end = "")
		elif x[n-i][0] < 0:
			if i-1 == 1:
				print( " " + str(round(x[n-i][0],3)) + "*x"  ,end = "")
			elif i-1 == 0:
				print( " " + str(round(x[n-i][0],3)))
			else:
				print( " " + str(round(x[n-i][0],3)) + "*x^" + str(i-1) ,end = "")
		elif x[n-i][0] == 1:
			if i-1 == 1:
				print( "x"  ,end = "")
			elif i-1 == 0:
				print( " +" + str(round(x[n-i][0],3)))
			elif i == n:
				print( "x^" + str(i-1) ,end = "")
			else:
				print(  " x^" + str(i-1) ,end = "")

a = [[1],[-12],[30],[0],[12]]
#x,f = initializeaza(a)
print("-----------------Interpolare Lagrange---------------")
interpolareLagrange(x_prim,x,f)

print("\n-----------metoda cele mai mici patrate-------------")

x,f = initializeaza(a)
print("  x :",x)
print("f(x):",f)
metoda(x_prim,x,f,4)

afiseazaPolinom(a)
print(d(len(a)-1,a,x_prim))
print(func(x_prim))

