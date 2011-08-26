#!/usr/bin/env python

from scipy import *

# Number of sequences distance hp from 0 and 1 from a sequence at h:
def gcond(hp,h,L):

	if hp == h-1:
		return h

	if hp == h:
		return 2*h

	if hp == h+1:
		return 3*(L-h)

	return 0
	

# Derivatives:
def ddt (Y,V, L, p):

	# Rate of mutation to neighbouring sequences:
	mup = p['mu']/(3*L)

	dYdt = (1.-p['mu'])*p['beta']*p['xbar']*V

	for h in range(len(Y)):

		dYdt[h] += mup*p['beta']*p['xbar']*gcond(h,h,L)*V[h]

		if h>0:
			dYdt[h] += mup*p['beta']*p['xbar']*gcond(h-1,h,L)*V[h-1]

		if h<L:
			dYdt[h] += mup*p['beta']*p['xbar']*gcond(h+1,h,L)*V[h+1]
	
	dVdt = p['k']*Y - p['beta']*p['xbar']*V

	return (dYdt,dVdt)

# Semi-implicit integration step:
def step(Y, V, dt, L, p):

	Yp = Y.copy()
	Vp = V.copy()

	for i in range(3):
		dYdt,dVdt = ddt(Y, V, L, p)
		Yp = Y + 0.5*dt*dYdt
		Vp = V + 0.5*dt*dVdt
	
	Y = 2.*Yp - Y
	V = 2.*Vp - V

	return (Y, V)

# Send sample to stdout:
def sample(V, t):
	print t,
	for i in range(len(V)):
		print V[i],
	print

# Produce header for sample file:
def header(V):
	print 't',
	for i in range(len(V)):
		print 'n'+str(i),
	print


## MAIN ##
if __name__ == '__main__':

	# Model parameters:
	p = {}
	L = 105
	p['beta'] = 5e-13
	p['k'] = 1e3
	p['xbar'] = 2.5e11
	p['mu'] = 2e-5*L

	# Simulation parameters:
	T = 2.
	Nt = 1001
	Nsamples = 101

	# Derived simulation parameters:
	dt = T/(Nt-1)
	steps_per_sample = (Nt-1)/(Nsamples-1)

	# Initial condition:

	Y = zeros(L+1)
	V = zeros(L+1)
	V[0] = 1

	# Generate header:
	header(V)

	# Perform first sample:
	sample(V, 0)

	# Integration loop:
	for tidx in range(Nt):

		# Perform semi-implicit integration step:
		Y,V = step(Y, V, dt, L, p)

		# Sample if necessary:
		if tidx % steps_per_sample == 0:
			sample(V, tidx*dt)
