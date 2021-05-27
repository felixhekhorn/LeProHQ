# LeProHQpy

(Un-)Polarized Lepto-Production of Heavy Quarks.

This is the stand-alone Python wrapper for the fully-inclusive coefficient functions.

To see this implementation of the coefficient functions in action, i.e. actual structure functions, please use [yadism](https://n3pdf.github.io/yadism/).

## Normalization

The available structure functions are
$$
F_2, F_L, xF_3, 2xg_1, g_4, g_L
$$
as defined by the [PDG](https://pdg.lbl.gov/).

The normalization of the factorization formula is given by
$$
F_2^Q(x,Q^2) = \frac{\alpha_s \xi}{4\pi^2} \sum\limits_{j=q,\overline{q},g} \int\limits_x^{z_{max}} \frac{dz}{z} f_j(x/z,\mu_F^2) c_{j}(\xi, \eta)
$$
with $\xi = Q^2/m^2$ and $\eta = (s-4m^2)/(4m^2)$ as scaling variables and $z_{max}=Q^2/(4m^2+Q^2)$ the kinematic bound for the final state.
The partonic coefficient functions are then given by
$$
c_g = e_Q^2 \left( c_g^{(0)} + 4\pi\alpha_s\left( c_g^{(1)} + \overline c_g^{(1),F} \ln(\mu_F^2/m^2) + \overline c_g^{(1),R} \ln(\mu_R^2/m^2) \right) \right)
$$
and
$$
c_q = 4\pi\alpha_s \left( e_Q^2\left( c_q^{(1)} + \overline c_q^{(1),F} \ln(\mu_F^2/m^2) \right) + e_q^2 d_q^{(1)} + \right)
$$
in the case of electroproduction. The extension to the case of leptoproduction can be obtained by acting accordingly on the charge factors.

