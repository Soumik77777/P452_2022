import numpy as np

n=3
def right(x) :
	if x%n==0 :
		return x-n-1
	else :
		return x+1

def left(x) :
	if (x-1)%n==0 :
		return x+n-1
	else :
		return x-1

def up(x) :
	if (x-1)//n < 1 :
		return n*(n-1)+x
	else :
		return x-n

def don(x) :
	if (x-1)//n>1 :
		return (((x+n)//n)-n)+((x+n)%n)
	else :
		return x+n


def Gen(i,j) :					## Function to generate A "on the fly"
	if i==j :
		return 0.96
	elif right(i+1)==j+1 :
		return 0.5
	elif left(i+1)==j+1 :
		return 0.5
	elif up(i+1)==j+1 :
		return 0.5
	elif don(i+1)==j+1 :
		return 0.5
	else :
		return 0

print (Gen(1,2))










############### Storing "if" needed




A=[[0 for i in range (n**2)] for j in range (n**2)]
print(len(A))
for i in range (n**2) :
	for j in range (n**2) :
		if i==j :
			A[i][j]=0.96
		elif right(i+1)==j+1 :
			A[i][j]=0.5
		elif left(i+1)==j+1 :
			A[i][j]=0.5
		elif up(i+1)==j+1 :
			A[i][j]=0.5
		elif don(i+1)==j+1 :
			A[i][j]=0.5
		else :
			A[i][j]=0

for i in range(len(A)):
	print(A[i])


'''
0.5
9
[0.96, 0.5, 0.5, 0.5, 0, 0, 0.5, 0, 0]
[0.5, 0.96, 0.5, 0, 0.5, 0, 0, 0.5, 0]
[0, 0.5, 0.96, 0, 0, 0.5, 0, 0, 0.5]
[0.5, 0, 0, 0.96, 0.5, 0.5, 0.5, 0, 0]
[0, 0.5, 0, 0.5, 0.96, 0.5, 0, 0.5, 0]
[0, 0.5, 0.5, 0, 0.5, 0.96, 0, 0, 0.5]
[0.5, 0, 0, 0.5, 0, 0, 0.96, 0.5, 0.5]
[0, 0.5, 0, 0, 0.5, 0, 0.5, 0.96, 0.5]
[0.5, 0, 0, 0, 0.5, 0.5, 0, 0.5, 0.96]
'''
