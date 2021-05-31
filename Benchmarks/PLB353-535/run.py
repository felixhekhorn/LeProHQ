import LeProHQpp

nlf = 3
m2 = 2.
Delta = 1e-6

obj = LeProHQpp.InclusiveLeptoProduction(nlf,m2,Delta)
obj.setQ2(10.)
obj.setPdf("DSSV14",0)
obj.setAlphaSByLHAPDF("MSTW2008nlo90cl",0)
obj.setXBjorken(.1)
obj.flags().useNextToLeadingOrder = False
print(obj.F())