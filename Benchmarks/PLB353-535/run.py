import LeProHQpp

nlf = 3
m2 = 2.0
Delta = 1e-6

obj = LeProHQpp.InclusiveLeptoProduction(nlf, m2, Delta)
obj.setQ2(10.0)
obj.setPdf("GRV94LO", 0)
obj.setAlphaSByLHAPDF("MSTW2008nlo90cl", 0)
obj.setXBjorken(0.1)
obj.flags().useNextToLeadingOrder = False
print(obj.F())
