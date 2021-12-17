#!/usr/bin/python
# -*- coding: utf-8 -*-
# fmt: off

import Util

import numpy as np

from ElProduction import projT, OrderFlag_LO, OrderFlag_NLO, ChannelFlag_Full
from InclusiveRunner import InclusiveRunner

path = Util.getAbsProjectPath()+"data/NPB540-345/"
Util.pInfo("Computing files for NPB540-345 ...")

r = InclusiveRunner();

m2 = 1.5*1.5
q2 = -1e-2
lambdaQCD = 0.2

Delta = 1e-6

# add data point
def addXtOfGRSV96STD(sqrtSh,n,cMuF2,cMuR2):
	pdf = "GRSV96STD"+("" if 0 == n else "N")+"LO";
	
	a = (pdf,sqrtSh,cMuF2*100.,cMuR2*100.,n)
	fpXt = path+"dFP_dxt-inc-%s-sqrtSh_%g-cMuF2_%g-cMuR2_%g-%d.dat"%a
	
	N = 100
	for j in xrange(N-1):
		xt = 10**(-3. + 3./N*(.5 + j))
		r.add({
			"objArgs":(m2,q2,Delta, projT.P,3,),
			"pdf": (pdf,0,),
			"lambdaQCD": lambdaQCD,
			"muF2": (1.*cMuF2,0.,0.,1.*cMuF2),
			"muR2": (1.*cMuR2,0.,0.,1.*cMuR2),
			"hadronicS":sqrtSh*sqrtSh,
			"f": ("dF_dHAQTransverseMomentumScaling", xt, (OrderFlag_LO if 0 == n else OrderFlag_NLO), ChannelFlag_Full),
			"msg": "dxt %s sqrtSh=%g, cMuF2=%g, cMuR2=%g, n=%g (xt = %e)"%(pdf,sqrtSh,cMuF2,cMuR2,n,xt),
			"fp": fpXt, "var": xt
		})
def addFig10():
	for n in [0]:#[0,1]:
		for cMuF2 in [.25,.5,1.,2.,4.]:
			for cMuR2 in [.25,.5,1.,2.,4.]:
				addXtOfGRSV96STD(10.,n,cMuF2,cMuR2)
	addXtOfGRSV96STD(200,0,1.,1.)
	#addXtOfGRSV96STD(200,1,1.,1.)

# add data point
def addFig11Helper(proj,pdf,n):
	sqrtSh = 10.
	a = (proj,pdf,sqrtSh,n)
	fpRap = path+"dF%s_dy-inc-%s-sqrtSh_%g-%d.dat"%a
	fpXf = path+"dF%s_dxF-inc-%s-sqrtSh_%g-%d.dat"%a
	N = 100
	for j in xrange(N-1):
		y = -2. + 4./N*(.5 + j)
		r.add({
			"objArgs":(m2,q2,Delta,proj,3,),
			"pdf": (pdf,0,),
			"lambdaQCD": lambdaQCD,
			"mu2": (2.,0.,0.,0.),
			"hadronicS":sqrtSh*sqrtSh,
			"f": ("dF_dHAQRapidity", y, (OrderFlag_LO if 0 == n else OrderFlag_NLO), ChannelFlag_Full),
			"msg": "dy %s sqrtSh=%g, n=%g (y = %e)"%(pdf,sqrtSh,n,y),
			"fp": fpRap, "var": y
		})
		xF = -1. + 2./N*(.5 + j)
		r.add({
			"objArgs":(m2,q2,Delta,proj,3,),
			"pdf": (pdf,0,),
			"lambdaQCD": lambdaQCD,
			"mu2": (2.,0.,0.,0.),
			"hadronicS":sqrtSh*sqrtSh,
			"f": ("dF_dHAQFeynmanX", xF, (OrderFlag_LO if 0 == n else OrderFlag_NLO), ChannelFlag_Full),
			"msg": "dxF %s sqrtSh=%g, n=%g (xF = %e)"%(pdf,sqrtSh,n,xF),
			"fp": fpXf, "var": xF
		})
def addFig11():
	addFig11Helper(projT.P,"GRSV96STDLO",0)
	#addFig11Helper(projT.P,"GRSV96STDNLO",1)
	#addFig11Helper(projT.G,"GRV94NLO",1)

# run
addFig10()
#addFig11()
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
