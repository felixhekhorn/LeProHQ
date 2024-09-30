import numpy as np
from scipy.interpolate import CubicSpline

# TODO: move this into a data file
logxis = np.array([-6.,-5.9,-5.8,-5.7,-5.6,-5.5,-5.4,-5.3,-5.2,-5.1,-5.,-4.9,-4.8,-4.7,-4.6,-4.5,-4.4,-4.3,-4.2,-4.1,-4.,-3.9,-3.8,-3.7,-3.6,-3.5,-3.4,-3.3,-3.2,-3.1,-3.,-2.9,-2.8,-2.7,-2.6,-2.5,-2.4,-2.3,-2.2,-2.1,-2.,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.])
vals = {("F2", "VV"): np.array([4.31775e-6, 5.36701e-6, 6.67017e-6, 8.28834e-6,
0.0000102973, 0.0000127909, 0.0000158856, 0.0000197252,
0.0000244882, 0.0000303952, 0.0000377195, 0.000046799,
0.0000580514, 0.0000719934, 0.0000892633, 0.00011065,
0.000137127, 0.000169897, 0.000210444, 0.000260599, 0.000322616,
0.000399279, 0.000494012, 0.000611035, 0.000755539, 0.000933909,
0.00115399, 0.00142544, 0.00176009, 0.00217248, 0.00268042,
0.00330576, 0.00407523, 0.00502155, 0.00618472, 0.00761361,
0.00936784, 0.0115201, 0.014159, 0.0173921, 0.0213505,
0.0261929, 0.032112, 0.0393408, 0.048161, 0.0589126, 0.072005,
0.0879306, 0.10728, 0.130762, 0.159221, 0.193666, 0.235296,
0.285534, 0.34606, 0.418859, 0.506261, 0.610996, 0.736251,
0.885731, 1.06373, 1.27518, 1.52578, 1.82201, 2.17126,
2.58189, 3.06332, 3.62608, 4.28195, 5.04395, 5.92647, 6.94529,
8.11761, 9.46211, 10.999, 12.7498, 14.7379, 16.9877, 19.5254,
22.3784, 25.5756, 29.1471, 33.1242, 37.5395, 42.4266, 47.8201,
53.7557, 60.2696, 67.3993, 75.1826, 83.6582, 92.8653, 102.844,
113.633, 125.275, 137.81, 151.279, 165.724, 181.187, 197.709,
215.334, 234.104, 254.061, 275.247, 297.706, 321.481, 346.615,
373.15, 401.13, 430.598, 461.596, 494.169, 528.36, 564.211,
601.767, 641.07, 682.163, 725.091, 769.897, 816.623, 865.313,
916.011, 968.76, 1023.6, 1080.58, 1139.75, 1201.13, 1264.79,
1330.75, 1399.07, 1469.79, 1542.95, 1618.6, 1696.77, 1777.52,
1860.88, 1946.89, 2035.62, 2127.08, 2221.34, 2318.42]),
            ("FL", "VV"): np.array([2.37037e-7, 2.98412e-7, 3.75678e-7, 4.7295e-7,
5.95409e-7, 7.49575e-7, 9.43658e-7, 1.18799e-6,
1.4956e-6, 1.88284e-6, 2.37035e-6, 2.98409e-6,
3.75674e-6, 4.72944e-6, 5.954e-6, 7.4956e-6,
9.43636e-6, 0.0000118796, 0.0000149554, 0.0000188276,
0.0000237022, 0.0000298389, 0.0000375642, 0.0000472895,
0.0000595323, 0.0000749442, 0.0000943452, 0.000118767,
0.00014951, 0.000188207, 0.000236916, 0.000298225, 0.00037539,
0.000472505, 0.000594722, 0.000748515, 0.000942025, 0.00118548,
0.00149173, 0.00187689, 0.00236122, 0.00297008, 0.00373526,
0.00469656, 0.00590371, 0.00741883, 0.00931929, 0.0117014,
0.0146845, 0.0184165, 0.0230797, 0.0288981, 0.0361457,
0.0451558, 0.0563314, 0.0701561, 0.0872054, 0.108158, 0.133802,
0.165048, 0.202921, 0.248565, 0.303222, 0.368217, 0.444916,
0.534693, 0.63887, 0.758666, 0.895142, 1.04915, 1.22128,
1.41187, 1.62094, 1.84828, 2.09337, 2.3555, 2.63378, 2.92716,
3.23451, 3.55466, 3.88642, 4.22861, 4.58009, 4.9398, 5.30674,
5.67998, 6.05871, 6.44217, 6.82969, 7.22069, 7.61466, 8.01114,
8.40975, 8.81016, 9.21208, 9.61526, 10.0195, 10.4246, 10.8305,
11.237, 11.6439, 12.0513, 12.4591, 12.8671, 13.2754, 13.6839,
14.0925, 14.5012, 14.9101, 15.3191, 15.7281, 16.1372, 16.5463,
16.9555, 17.3647, 17.7739, 18.1832, 18.5924, 19.0017, 19.411,
19.8203, 20.2297, 20.639, 21.0483, 21.4576, 21.867, 22.2763,
22.6856, 23.095, 23.5043, 23.9137, 24.323, 24.7324, 25.1417,
25.5511, 25.9604, 26.3697, 26.7791, 27.1884, 27.5978, 28.0071]),
            ("x2g1","VV"): np.array([3.55556e-7, 4.47618e-7, 5.63518e-7, 7.09427e-7,
8.93115e-7, 1.12437e-6, 1.41549e-6, 1.782e-6,
2.2434e-6, 2.82428e-6, 3.55555e-6, 4.47618e-6,
5.63517e-6, 7.09426e-6, 8.93115e-6, 0.0000112436,
0.0000141549, 0.00001782, 0.000022434, 0.0000282427,
0.0000355555, 0.0000447616, 0.0000563515, 0.0000709423,
0.0000893109, 0.000112436, 0.000141548, 0.000178198,
0.000224337, 0.000282422, 0.000355546, 0.000447603, 0.000563494,
0.000709389, 0.000893055, 0.00112427, 0.00141534, 0.00178176,
0.00224303, 0.00282368, 0.00355461, 0.00447468, 0.00563279,
0.00709049, 0.00892518, 0.0112342, 0.01414, 0.0177963,
0.0223966, 0.0281837, 0.0354622, 0.0446144, 0.0561193,
0.0705763, 0.0887348, 0.11153, 0.140126, 0.17597, 0.220853,
0.276986, 0.347087, 0.43448, 0.543207, 0.678151, 0.845171,
1.05123, 1.30455, 1.61473, 1.99284, 2.45159, 3.00533, 3.67016,
4.4639, 5.40617, 6.51828, 7.82323, 9.34566, 11.1118, 13.1492,
15.487, 18.1556, 21.1866, 24.6128, 28.4679, 32.7867, 37.6051,
42.9594, 48.887, 55.4258, 62.6146, 70.4924, 79.099, 88.4746,
98.6596, 109.695, 121.622, 134.483, 148.318, 163.171, 179.082,
196.096, 214.253, 233.598, 254.171, 276.018, 299.18, 323.7,
349.622, 376.988, 405.842, 436.227, 468.186, 501.763, 537.001,
573.942, 612.631, 653.111, 695.425, 739.616, 785.728, 833.805,
883.889, 936.024, 990.253, 1046.62, 1105.17, 1165.94, 1228.98,
1294.33, 1362.04, 1432.14, 1504.69, 1579.72, 1657.28, 1737.41,
1820.16, 1905.56, 1993.67, 2084.52, 2178.16, 2274.64])
}

for proj in ("F2", "FL","x2g1"):
    vals[(proj,"AA")] = vals[(proj, "VV")]
vals[("xF3","VA")] =vals[("x2g1","VV")]
vals[("g4","VA")] = -vals[("F2","VV")]
vals[("gL","VA")] = -vals[("FL","VV")]

        
def Adler(proj, cc, xi):
    """Adler sum rule correction term.
    
    This is the integral over the NLO Coulomb quark scaling coefficient functions."""
    return CubicSpline(logxis, vals[(proj, cc)])([np.log10(xi)])[0]