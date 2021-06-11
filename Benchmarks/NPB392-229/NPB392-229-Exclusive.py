#!/usr/bin/python
# -*- coding: utf-8 -*-

import Util

import numpy as np

from ElProduction import projT, ExclusiveHistT, OrderFlag_LO, OrderFlag_NLO, ChannelFlag_Full
from ExclusiveRunner import ExclusiveRunner

#pdf = "MSTW2008nlo90cl"
pdf = "MorfinTungB"
pathOut = Util.getAbsProjectPath()+"data/NPB392-229-%s/"%pdf
Util.pInfo("Computing files for NPB392-229 with "+pdf+" ...")

r = ExclusiveRunner();

xTilde = .8
omega = 1.
deltax = 1e-6
deltay = 7e-6
q2 = -10.

# add data point
def add(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,activatedHistograms,n,msg):
	its = 5 if 0 == n else 1
	calls = 500000 if 0 == n else 500000
	r.add({
		"objArgs":(m2,q2,proj,nlf,xTilde,omega,deltax,deltay,),
		"pdf": (pdf,0,), "lambdaQCD": lambdaQCD, "mu2": mu2, "bjorkenX":bjorkenX, 
		"activatedHistograms": activatedHistograms, "orderFlag": (OrderFlag_LO if 0 == n else OrderFlag_NLO), "channelFlag": ChannelFlag_Full,
		"calls": calls, "iterations": its, "verbosity": 2, "adaptChi2": 0 == n,
		"msg": msg
	})
def addPt(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,ptmax,n):
	qL = "c" if 3 == nlf else ("b" if 4 == nlf else "t")
	fpPt = pathOut+"dF%s%s_dpt_x-%g_%d.dat"%(proj,qL,-np.log10(bjorkenX),n)
	activatedHistograms = [
		(ExclusiveHistT.HAQTransverseMomentum,100,fpPt,0.,ptmax,)
	]
	add(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,activatedHistograms,n,fpPt)
def addRap(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,y0,n):
	qL = "c" if 3 == nlf else ("b" if 4 == nlf else "t")
	fpRap = pathOut+"dF%s%s_dy_x-%g_%d.dat"%(proj,qL,-np.log10(bjorkenX),n)
	activatedHistograms = [
		(ExclusiveHistT.HAQRapidity,100,fpRap,-y0,y0,)
	]
	add(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,activatedHistograms,n,fpRap)
# add charm data points
def addPtCharm(proj,bjorkenX,ptmax):
	m2 = 1.5**2
	nlf = 3
	lambdaQCD = 0.194
	mu2 = (4.,-1.,0.,4.,)
	addPt(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,ptmax,0)
	#addPt(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,ptmax,1)
def addRapCharm(proj,bjorkenX,y0):
	m2 = 1.5**2
	nlf = 3
	lambdaQCD = 0.194
	mu2 = (4.,-1.,0.,0.,)
	addRap(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,y0,0)
	#addRap(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,y0,1)
# add plots
def addPtCharmG():
	addPtCharm(projT.G, .1,     5.)
	addPtCharm(projT.G, .01,    8.)
	addPtCharm(projT.G, .001,  15.)
	addPtCharm(projT.G, .0001, 20.)
def addPtCharmL():
	addPtCharm(projT.L, .1,     5.)
	addPtCharm(projT.L, .01,    8.)
	addPtCharm(projT.L, .001,  15.)
	addPtCharm(projT.L, .0001, 20.)
def addRapCharmG():
	addRapCharm(projT.G, .1,    2. )
	addRapCharm(projT.G, .01,   3.5)
	addRapCharm(projT.G, .001,  4.5)
	addRapCharm(projT.G, .0001, 5.5)
def addRapCharmL():
	addRapCharm(projT.L, .1,    2. )
	addRapCharm(projT.L, .01,   3.5)
	addRapCharm(projT.L, .001,  4.5)
	addRapCharm(projT.L, .0001, 5.5)

addPtCharmG()
addPtCharmL()
addRapCharmG()
addRapCharmL()
r.run()

