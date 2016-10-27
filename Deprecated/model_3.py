#!/usr/bin/env python 

'''
                    Space Telescope Science Institute

Synopsis:  

Try to create a model for persistence in a diode


Command line usage (if any):

	usage: model.py filename

Description:  

Primary routines:

Notes:
									   
History:

120424 ksl Coding begun

'''

import sys
import math
import numpy
import pylab


def gen_ramps(rate,tint,nint,ttot):
	'''
	Generate a time history for the exposre 
	as a fundtion of time
	'''
	if ttot<tint*nint:
		ttot=tint*nint
	
	t=numpy.arange(0,ttot,1)
	x=numpy.zeros(ttot)
	i=0
	while i<nint*tint:
		z=i % tint
		print('z',z, rate)
		x[i]=z*rate
		i=i+1
	return t,x

def init_diode():
	'''
	Initialize the trap distribution.  Assume
	a gaussian distribution of traps centered
	somewhere around saturation

	Begin with none of the traps filled
	'''
	from scipy.stats import norm

	nlevels=100
	nmax=2e5
	saturation=80000
	width=10000
	ntraps=1e6
	nmin=nmax/nlevels
	# levels=numpy.arange(nmin,nmax,nlevels)
	levels=numpy.linspace(nmin,nmax,nlevels)
	density=norm.pdf(levels,saturation,width)
	density=density*ntraps
	frac=numpy.zeros_like(density)
	# print len(levels),len(density)
	return levels,density,frac


def pdf_pow(t,tmin,tmax,gamma):
	'''
	Calculate the normalized pdf for a power
	law of the type
		prob(t) propto t**-gamma
	between tmin and tmax allowing for the possibility 
	that gamma is 1
	'''

	if t<tmin:
		return 0
	if t>tmax:
		return 0
	if gamma==1:
		c=1./math.log(tmax/tmin)
		return c/t
	else:
		c=(-gamma+1)/(tmax**(-gamma+1)-tmin**(-gamma+1))
		return c*t**-gamma

def test_pdf_pow(gamma=2,tmin=10.,tmax=100.):
	'''
	Test the generation power density function of a power law
	'''


	x=numpy.linspace(tmin,tmax,tmax-tmin)

	sum=0
	for one in x:
		sum=sum+pdf_pow(one,10,100,gamma)

	print(sum)



def xinit_diode():
	'''
	Initialize the trap distribution.  Assume
	a gaussian distribution of traps centered
	somewhere around saturation

	Begin with none of the traps filled

	Here	levels are effectively 'energy' levels
		times  are the release times. 
	
	The dislocation are assumed to go accourding to
	a power law
	'''
	from scipy.stats import norm

	nlevels=100
	nmax=2e5
	saturation=80000
	width=10000
	ntraps=1e6
	nmin=nmax/nlevels
	# levels=numpy.arange(nmin,nmax,nlevels)
	levels=numpy.linspace(nmin,nmax,nlevels)

	ntimes=10
	tmin=100.
	tmax=10000.
	gamma=1.
	tconstants=numpy.linspace(tmin,tmax,ntimes)
	tstep=tconstants[1]-tconstants[0]
	# tconstants is an array containg the time constants

	print('test',len(tconstants))
	print('tconstants',tconstants)


	density=norm.pdf(levels,saturation,width)
	density=density*ntraps
	# density is one-d arraty containing the total number of traps at a level

	tdist=[]
	for one in tconstants:
		tdist.append(pdf_pow(one,tmin,tmax,gamma)*tstep)
	print('Time distribtion',tdist)
	tdist=numpy.array(tdist)
	# tdist is the probability that a trap will have a given time constant

	i=0
	while i < len(density):
		xxx=density[i]*tdist
		if i==0:
			density2=xxx
		else:
			density2=numpy.row_stack((density2,xxx))
		i=i+1
	

	frac=numpy.zeros_like(density2)
	return levels,tconstants,density2,frac

def dndt(x,levels,frac,up=0.100,down=0.02):
	'''
	Change the fraction of the levels that
	are occupied accrding to a simple model
	'''
	i=0
	z=[]
	while i<len(levels):
		if x<levels[i]:
			# change=-down*frac[i]
			change=-down*frac[i]**2.
			# change=-down*frac[i]**6
			z.append(change)
		else:
			change=up*(1.-frac[i]) 
			z.append(change)
		i=i+1
	z=numpy.array(z)
	z=frac+z
	return z



def xdndt(x,levels,tconstants,frac,ratio):
	'''
	Change the fraction of the levels that
	are occupied accrding to a simple model

	Note that is the 2d model where there
	are multiple time constants.

	frac is the fraction of the level population
	that is filed at each time constand

	The program returns the fraction of 
	traps that are populated in a second
	at the electron level x

	120525	ksl	Coded an verified that it 
			seems to preform as one
			would expect
	
	'''
	i=0
	z=[]
	while i<len(levels):
		if x<levels[i]:
			# print 'down'
			change=-frac[i]/tconstants
			# print change
		else:
			# print 'up'
			change=(1.-frac[i])/(ratio*tconstants)
			print(change)
		z.append(change)
		i=i+1
	z=numpy.array(z)
	z=frac+z
	return z


def plot_frac(levels,frac):
	'''
	'''
	pylab.figure(8)

	y=[]
	for one in frac:
		y.append(numpy.average(one))
	pylab.plot(levels,y,'-')
	pylab.draw()
	return




def current(tstart,t,tot):
	'''
	Plot the current given a starting time
	in log log coordinates
	'''
	z=t-tstart
	i=0
	x=[]
	y=[]
	while i < len(z):
		if z[i]>0:
			x.append(z[i])
			y.append(tot[i]-tot[i-1])
		i=i+1
	
	x=numpy.array(x,dtype=float)
	y=numpy.array(y,dtype=float)
	y=-y
	# print x
	# print y

	mid=len(x)/2
	q=x/x[mid]

	# print 'test did red'

	# print q

	

	z=y[mid]*(x/x[mid])**-1
	# print z

	pylab.figure(5,(6,6))
	pylab.clf()
	pylab.loglog(x,y,'o')
	pylab.loglog(x,z,'r-')
	# print 'test did red'

	
	pylab.xlabel('Time(s)')
	pylab.ylabel('Persistence (e/s)')
	pylab.title('Persistence')
	pylab.draw()
	pylab.savefig('persistence_all.png')
	z=pylab.xlim()
	pylab.xlim(100,z[1])
	pylab.axis((100,10000,0.01,1))
	pylab.draw()
	pylab.savefig('persistence.png')
		
def do_spec(time=1000,tint=100,nint=3,up=0.03,down=0.03):
	'''
	Calculate the persistence spectrum at a specific
	point in time
	'''
	# Generate a ramp assuming a rate of 1/tint so at the end of
	# the ramp the total is 1
	# Here we end the ramp at time we want to read everything out
	t,x=gen_ramps(1./tint,tint,nint,time+tint*nint)
	levels,density,frac=init_diode()

	print(x)

	# Generate a numpy array from 10,000 electrons to 1e6 electons
	# Do it logarithmically
	values=numpy.logspace(4,6,21)

	result=[]
	for one in values:
		xx=x*one
		xlevels=levels
		xfrac=frac
		i=0
		while i<len(t):
			# Calculate the fraction of traps that are
			# filed as a function of time
			# if i % 10 == 0:
			# 	print 'test',xx[i]
			xfrac=dndt(xx[i],xlevels,xfrac,up,down)
			# calulate the number of filled traps 
			zz=density*xfrac
			tot=numpy.average(zz)*len(zz)
			# tot heare is the number of filled traps
			i=i+1
		result.append(tot)
		print('ok',one,tot)
	print(result)







def xdoit(rate=1000,tint=100,nint=3,ttot=10000,ratio=10.):
	'''
	Mimic a ramp and the persistence that results, where

	rate is the rate at which photons are detected during a ramp
	tint is the length of a ramp in seconds
	nint is the number of repreates
	ttot is the total amount of time to be simulated
	ratio is the ratio of trapping vs release times

	120522 These particular choices give a pretty good power law decay between 100 and 10**5 
	seconds, given the current model for dn/dt
	120525	This is the version of the model that attempts to allow for varying trap times
	'''

	# Generate an exposure history
	t,x=gen_ramps(rate,tint,nint,ttot)
	pylab.figure(1,(6,6))
	pylab.clf()
	pylab.plot(t,x,'-')
	pylab.xlabel('Time(s)')
	pylab.ylabel('Well depth (e)')
	pylab.title('Ramp history')
	pylab.xlim(0,tint*(nint+2))
	pylab.draw()
	pylab.savefig('exposure_hist.png')


	# Generate the trap distribution and start with
	# all of the traps empty
	# levels,density,frac=init_diode()
	levels,tconstants,density2,frac=xinit_diode()

	# For the 2d modeld need to sum the density to fix figure1

	density=[]
	for one in density2:
		sum=one.sum()
		density.append(sum)
	density=numpy.array(density)


	pylab.figure(2,(6,6))
	pylab.clf()
	pylab.plot(levels,density,'-')
	pylab.xlabel('Well depth')
	pylab.ylabel('Available traps')
	pylab.title('Trap distibution')
	pylab.draw()
	pylab.savefig('trap_distrib.png')

	pylab.figure(3,(6,6))
	pylab.clf()

	i=0
	xtot=[]
	while i<len(t):
		# Calculate the fraction of traps that are
		# filed as a function of time
		# frac=xdndt(x[i],levels,frac,up,down)
		frac=xdndt(x[i],levels,tconstants,frac,ratio)
		
		# calulate the number of filled traps 
		zz=density2*frac
		tot=numpy.average(zz)*len(zz)
		# tot heare is the number of filled traps
		xtot.append(tot)
		if t[i] % 50 == 0:
			z=[]
			for one in frac:
				z.append(one.sum())
			pylab.plot(levels,z,'-')
			pylab.draw()
			print(t[i],tot)
		i=i+1


	pylab.title('Level fractions')
	pylab.xlabel('Well depth')
	pylab.draw()
	pylab.savefig('level_frac.png')

	pylab.figure(4,(6,6))
	pylab.clf()
	pylab.plot(t,xtot,'-')
	pylab.xlabel('Time(s)')
	pylab.ylabel('Filled traps')
	pylab.title('Filled traps with time')
	pylab.xlim(0,tint*(nint+2))
	pylab.draw()
	pylab.savefig('traps_vs_time.png')

	current(tint*nint,t,xtot)


	



def doit(rate=1000,tint=100,nint=3,ttot=10000,up=0.03,down=0.03):
	'''
	Mimic a ramp and the persistence that results, where

	rate is the rate at which photons are detected during a ramp
	tint is the length of a ramp in seconds
	nint is the number of repreates
	ttot is the total amount of time to be simulated
	up   is the rate at which holes are filed when expposed
	down is the normalization factor for the rate decay

	120522 These particular choices give a pretty good power law decay between 100 and 10**5 
	seconds, given the current model for dn/dt
	'''

	# Generate an exposure history
	t,x=gen_ramps(rate,tint,nint,ttot)
	pylab.figure(1,(6,6))
	pylab.clf()
	pylab.plot(t,x,'-')
	pylab.xlabel('Time(s)')
	pylab.ylabel('Well depth (e)')
	pylab.title('Ramp history')
	pylab.xlim(0,tint*(nint+2))
	pylab.draw()
	pylab.savefig('exposure_hist.png')


	# Generate the trap distribution and start with
	# all of the traps empty
	levels,density,frac=init_diode()

	pylab.figure(2,(6,6))
	pylab.clf()
	pylab.plot(levels,density,'-')
	pylab.xlabel('Well depth')
	pylab.ylabel('Available traps')
	pylab.title('Trap distibution')
	pylab.draw()
	pylab.savefig('trap_distrib.png')

	pylab.figure(3,(6,6))
	pylab.clf()

	i=0
	xtot=[]
	while i<len(t):
		# Calculate the fraction of traps that are
		# filed as a function of time
		frac=dndt(x[i],levels,frac,up,down)
		# calulate the number of filled traps 
		zz=density*frac
		tot=numpy.average(zz)*len(zz)
		# tot heare is the number of filled traps
		xtot.append(tot)
		if t[i] % 10 == 0:
			pylab.plot(levels,frac,'-')
			print(t[i],tot)
		i=i+1


	pylab.plot(levels,frac,'-')
	pylab.title('Level fractions')
	pylab.xlabel('Well depth')
	pylab.draw()
	pylab.savefig('level_frac.png')

	pylab.figure(4,(6,6))
	pylab.clf()
	pylab.plot(t,xtot,'-')
	pylab.xlabel('Time(s)')
	pylab.ylabel('Filled traps')
	pylab.title('Filled traps with time')
	pylab.xlim(0,tint*(nint+2))
	pylab.draw()
	pylab.savefig('traps_vs_time.png')

	current(tint*nint,t,xtot)




# Next lines permit one to run the routine from the command line
if __name__ == "__main__":
	import sys
	if len(sys.argv)>1:
		# doit(int(sys.argv[1]))
		doit(sys.argv[1])
	else:
		print('usage: model.py filename')
