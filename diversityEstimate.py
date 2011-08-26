#!/usr/bin/env python

from scipy import *

# Number of sequences distance hp from 0 and 1 from a sequence at h:
def g(hp,h,L):

	if hp == h-1:
		return h

	if hp == h:
		return 2*h

	if hp == h+1:
		return 3*(L-h)

	return 0
	

# Derivatives:
def ddt (X, Y, V, L, p):

	# Rate of mutation to neighbouring sequences:
	mup = p['mu']/(3.*L)

	dXdt = p['lambda'] - p['d']*X
	for h in range(L+1):
		dXdt -= p['beta']*X*V[h]
	
	dYdt = (1.-p['mu'])*p['beta']*X*V - p['a']*Y
	for h in range(L+1):

		dYdt[h] += mup*p['beta']*X*g(h,h,L)*V[h]

		if h>0:
			dYdt[h] += mup*p['beta']*X*g(h-1,h,L)*V[h-1]

		if h<0:
			dYdt[h] += mup*p['beta']*X*g(h+1,h,L)*V[h+1]
	
	dVdt = p['k']*Y - p['beta']*X*V - p['u']*V

	return (dXdt, dYdt,dVdt)

# Semi-implicit integration step:
def step(X, Y, V, dt, L, p):

	Xp = X
	Yp = Y.copy()
	Vp = V.copy()

	for i in range(3):
		dXdt,dYdt,dVdt = ddt(X, Y, V, L, p)
		Xp = X + 0.5*dt*dXdt
		Yp = Y + 0.5*dt*dYdt
		Vp = V + 0.5*dt*dVdt
	
	X = 2.*Xp - X
	Y = 2.*Yp - Y
	V = 2.*Vp - V

	return (X, Y, V)

# Send sample to stdout:
def sample(X, Y, V, t):
	print t, X,
	for h in range(len(V)):
		print Y[h], V[h],
	print

# Produce header for sample file:
def header(L):
	print 't X',
	for h in range(L+1):
		print 'Y'+str(h), 'V'+str(h),
	print


## MAIN ##
if __name__ == '__main__':

	# Model parameters:
	p = {}
	L = 105
	p['lambda'] = 2.5e8
	p['k'] = 1e3
	p['beta'] = 5e-13
	p['mu'] = 2e-5*L
	p['d'] = 1e-3
	p['a'] = 1.
	p['u'] = 3.

	# Simulation parameters:
	T = 10.
	Nt = 1001
	Nsamples = 101

	# Derived simulation parameters:
	dt = T/(Nt-1)
	steps_per_sample = (Nt-1)/(Nsamples-1)

	# Initial condition:

	X = p['lambda']/p['d']  # Uninfected steady state
	Y = zeros(L+1)
	V = zeros(L+1)
	V[0] = 1

	# Generate header:
	header(L)

	# Perform first sample:
	sample(X, Y, V, 0)

	# Integration loop:
	for tidx in range(Nt):

		# Perform semi-implicit integration step:
		X,Y,V = step(X, Y, V, dt, L, p)

		# Sample if necessary:
		if tidx % steps_per_sample == 0:
			sample(X, Y, V, tidx*dt)
