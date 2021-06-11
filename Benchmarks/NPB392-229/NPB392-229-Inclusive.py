#!/usr/bin/python
# -*- coding: utf-8 -*-

import Util

import numpy as np

from ElProduction import projT, OrderFlag_LO, OrderFlag_NLO, ChannelFlag_Full
from InclusiveRunner import InclusiveRunner

#pdf = "MSTW2008nlo90cl"
pdf = "MorfinTungB"
pathOut = Util.getAbsProjectPath()+"data/NPB392-229-%s/"%pdf
Util.pInfo("Computing files for NPB392-229 with "+pdf+" ...")

r = InclusiveRunner();

Delta = 1e-7
q2 = -10
fpPtTemplate = "dF%s%s_dpt_inc_x-%g_%d.dat"
fpRapTemplate = "dF%s%s_dy_inc_x-%g_%d.dat"
lsTo2 = []

# add data point
def addPt(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,ptmax,n):
	qL = "c" if 3 == nlf else ("b" if 4 == nlf else "t")
	fpPt = pathOut+fpPtTemplate%("%s",qL,-np.log10(bjorkenX),n)
	if projT.G == proj:
		lsTo2.append(fpPt)
	fpPt = fpPt%proj
	N = 100
	for j in xrange(N-1):
		pt = ptmax/N*(.5 + j)
		r.add({
			"objArgs":(m2,q2,Delta,proj,nlf,),
			"pdf": (pdf,0,), "lambdaQCD": lambdaQCD, "mu2": mu2, "bjorkenX":bjorkenX,
			"f": ("dF_dHAQTransverseMomentum", pt, (OrderFlag_LO if 0 == n else OrderFlag_NLO), ChannelFlag_Full),
			"msg": fpPt+" (pt = %e)"%(pt),
			"fp": fpPt,
			"var": pt
		})
def addRap(m2,proj,nlf,lambdaQCD,mu2,bjorkenX,y0,n):
	qL = "c" if 3 == nlf else ("b" if 4 == nlf else "t")
	fpRap = pathOut+fpRapTemplate%("%s",qL,-np.log10(bjorkenX),n)
	if projT.G == proj:
		lsTo2.append(fpRap)
	fpRap = fpRap%proj
	N = 100
	for j in xrange(N-1):
		y = y0*(-1. + 2./N*(.5 + j))
		r.add({
			"objArgs":(m2,q2,Delta,proj,nlf,),
			"pdf": (pdf,0,), "lambdaQCD": lambdaQCD, "mu2": mu2, "bjorkenX":bjorkenX,
			"f": ("dF_dHAQRapidity", y, (OrderFlag_LO if 0 == n else OrderFlag_NLO), ChannelFlag_Full),
			"msg": fpRap+" (y = %e)"%(y),
			"fp": fpRap,
			"var": y
		})
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
# run
addPtCharmG()
addPtCharmL()
addRapCharmG()
addRapCharmL()
l = r.run()

# reorder
m = {}
for e in l:
	if not m.has_key(e["fp"]): m[e["fp"]] = []
	m[e["fp"]].append((e["var"],e["res"]))
# write
def writeDataToFile(fp,d):
	with open(fp,"w") as o:
		for l in d:
			o.write("% e\t% e\n"%l)
	Util.pSucc(fp)
for fp in m:
	d = m[fp]
	sd = sorted(d, key = lambda a: a[0])
	m[fp] = sd
	writeDataToFile(fp,sd)
# combine 2
for fpT in lsTo2:
	dG = m[fpT%"G"]
	dL = m[fpT%"L"]
	d2 = []
	for j in xrange(len(dG)):
		d2.append((dG[j][0],dG[j][1] + 1.5*dL[j][1]))
	writeDataToFile(fpT%"2",d2)
	
