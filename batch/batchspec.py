#!/usr/bin/python
from Numeric import *
import sys
import os

#spx = arange(0.1, 20, 0.1) 
spx = arange(0.75, 1.25, 0.005) 
spx = arange(0.7, 1.31, 0.01) 
spx = arange(8.6, 9.2 + 0.005, 0.005) 

N = int(sys.argv[1])
n = int(sys.argv[2])
pid = int(sys.argv[3])
progname = sys.argv[4].lstrip().rstrip()

print "Running batch", n, "out of", N, "using binary :" + progname + ":" 

prog = "time ./" + progname
preargs  = " "
postargs = " 0.75 0.005 1.25"
postargs = " 11.585 0.0001 11.605"

# split sp array in N pieces.
def split_array(spx, N):
	sa = []
	for i in range(0, N):
		a = []
		sa.append(a)
	k = 0 
	for i in range(0, len(spx)):
		sa[k].append(spx[i])
		k = int((k+1)%N)
	return sa

# execute run
sa = split_array(spx, N)

for p1 in sa[int(n-1)]:
	print "Running", p1, "from", sa[int(n-1)]
	exec_str = str(prog)+str(preargs)+str(p1)+str(postargs)
	print "Executing: ", exec_str
	os.system(exec_str)
	#os.system(str(prog)+str(preargs)+str(p1)+" "+str(p2)+" "+str(postargs))
