from LeProHQpp import Projection as proj

pt = (
    dict(num=1, proj_=proj.F2, x=0.1, ptmax=5.0, yrange=(1e-5, 1e-2)),
    dict(num=2, proj_=proj.F2, x=0.01, ptmax=8.0, yrange=(1e-4, 1e-1)),
    dict(num=3, proj_=proj.F2, x=0.001, ptmax=15.0, yrange=(1e-4, 1e-1)),
    dict(num=4, proj_=proj.F2, x=0.0001, ptmax=20.0, yrange=(1e-4, 1e-1)),
    dict(num=6, proj_=proj.FL, x=0.1, ptmax=5.0, yrange=(1e-6, 1e-3)),
    dict(num=7, proj_=proj.FL, x=0.01, ptmax=10.0, yrange=(1e-5, 1e-2)),
    dict(num=8, proj_=proj.FL, x=0.001, ptmax=15.0, yrange=(1e-5, 1e-2)),
    dict(num=9, proj_=proj.FL, x=0.0001, ptmax=20.0, yrange=(1e-5, 1e-2)),
)
"""list of configurations for each transverse momentum figure"""

y = (
    dict(num=11, proj_=proj.F2, x=0.1, y0=2.0, yrange=(1e-6, 1e-2)),
    dict(num=12, proj_=proj.F2, x=0.01, y0=3.5, yrange=(1e-5, 1e-1)),
    dict(num=13, proj_=proj.F2, x=0.001, y0=4.5, yrange=(1e-5, 1e-1)),
    dict(num=14, proj_=proj.F2, x=0.0001, y0=5.5, yrange=(1e-5, 1e-1)),
    dict(num=16, proj_=proj.FL, x=0.1, y0=2.0, yrange=(1e-7, 1e-3)),
    dict(num=17, proj_=proj.FL, x=0.01, y0=3.0, yrange=(1e-6, 1e-2)),
    dict(num=18, proj_=proj.FL, x=0.001, y0=4.5, yrange=(1e-6, 1e-2)),
    dict(num=19, proj_=proj.FL, x=0.0001, y0=5.5, yrange=(1e-5, 1e-1)),
)
"""list of configurations for each rapidity figure"""
