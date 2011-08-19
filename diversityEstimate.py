#!/usr/bin/env python

from scipy import *

# Number of sequences distance h from 0:
def g(h, L):
	return float((3.**h)*comb(L, h))

# Derivative of n:
def dndt (n, L, gam, mu):

	res = gam*n

	for h in range(1,len(n)):
		res[h] += mu*g(h, L)/g(h-1, L)*n[h-1]

	return res

# Semi-implicit integration step:
def step(n, dt, L, gam, mu):

	np = n.copy()

	for i in range(3):
		np = n + 0.5*dt*dndt(n, L, gam, mu)
	
	n = 2.*np - n

	return n

# Send sample to stdout:
def sample(n, t):
	print t,
	for i in range(len(n)):
		print n[i],
	print

# Produce header for sample file:
def header(n):
	print 't',
	for i in range(len(n)):
		print 'n'+str(i),
	print


## MAIN ##
if __name__ == '__main__':

	# Model parameters:
	L = 105
	beta = 5e-13
	k = 1e3
	x = 2.5e11
	gam = sqrt(beta*k*x)
	mu = 2e-5*L
	n0 = 1.

	# Simulation parameters:
	T = 2.
	Nt = 1001
	Nsamples = 101

	# Derived simulation parameters:
	dt = T/(Nt-1)
	steps_per_sample = (Nt-1)/(Nsamples-1)

	# Initialise state vector:
	n = zeros(L+1)
	n[0] = n0

	# Generate header:
	header(n)

	# Perform first sample:
	sample(n, 0)

	# Integration loop:
	for tidx in range(Nt):

		# Perform semi-implicit integration step:
		n = step(n, dt, L, gam, mu)

		# Sample if necessary:
		if tidx % steps_per_sample == 0:
			sample(n, tidx*dt)
