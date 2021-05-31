# configure
conf_data = configuration_data()
conf_data.set('VERSION_STR', meson.project_version())
# pass PDF support
if get_option('morfintungb')
    conf_data.set('HAVE_MORFINTUNGB', '1')
endif
if get_option('cteq3')
    conf_data.set('HAVE_CTEQ3', '1')
endif
if get_option('grv')
    conf_data.set('HAVE_GRV', '1')
endif
if get_option('grsv')
    conf_data.set('HAVE_GRSV', '1')
endif
if get_option('dssv')
    conf_data.set('HAVE_DSSV', '1')
endif
# write
configure_file(output : 'config-buildtime.hpp',
               configuration : conf_data)

configuration_inc = include_directories('.')

# add sources
sources = []

sources += [
    'Common/ME/BQED.cpp',
    'Common/AbstractFixedOrderLeptoProduction.cpp',
    'Common/AbstractIntKer.cpp',
    'Common/AbstractLeptoProduction.cpp',
    'Common/EnvUtils.cpp',
]

sources += [
    'FullyDiffLeptoProduction.cpp',
    'FullyDiff/ME/A1.cpp',
    'FullyDiff/ME/A1Counter.cpp',
    'FullyDiff/ME/A2.cpp',
    'FullyDiff/ME/A3.cpp',
    'FullyDiff/ME/R.cpp',
    'FullyDiff/ME/RCounterX.cpp',
    'FullyDiff/ME/RCounterXY.cpp',
    'FullyDiff/ME/RCounterY.cpp',
    'FullyDiff/ME/ROK.cpp',
    'FullyDiff/ME/RQED.cpp',
    'FullyDiff/ME/SV.cpp',
    'FullyDiff/IntKer.cpp',
    'FullyDiff/PhasespacePoint.cpp',
]

sources += [
    'FullyInclusiveLeptoProduction.cpp',
    'FullyInclusive/coeffs/cg0.cpp',
    'FullyInclusive/coeffs/cg0t.cpp',
    'FullyInclusive/coeffs/cg1_a10.cpp',
    'FullyInclusive/coeffs/cg1.cpp',
    'FullyInclusive/coeffs/cgBar1.cpp',
    'FullyInclusive/coeffs/cq1.cpp',
    'FullyInclusive/coeffs/cqBarF1.cpp',
    'FullyInclusive/coeffs/dq1.cpp',
    'FullyInclusive/IntKer.cpp',
]

sources += [
    'FullyInclusiveSoftResummedLeptoProduction.cpp',
    'FullyInclusiveSoftResummed/cg0t.cpp',
    'FullyInclusiveSoftResummed/Interpolation.cpp',
    'FullyInclusiveSoftResummed/IntKer.cpp',
    'FullyInclusiveSoftResummed/MellinFuncs.cpp',
    'FullyInclusiveSoftResummed/ResExp.cpp',
]

sources += [
    'InclusiveLeptoProduction.cpp',
    'Inclusive/ME/IntA1.cpp',
    'Inclusive/ME/IntA2.cpp',
    'Inclusive/ME/IntA3.cpp',
    'Inclusive/ME/IntR.cpp',
    'Inclusive/ME/SV.cpp',
    'Inclusive/IntKer.cpp',
]

sources += [
    'Pdf/PdfWrapper.cpp'
]
# add PDF source files
if get_option('morfintungb')
    sources += ['Pdf/mortungb.f']
endif
if get_option('cteq3')
    sources += ['Pdf/cteq3.f']
endif
if get_option('grv')
    sources += ['Pdf/grv94.f']
endif
if get_option('grsv')
    sources += ['Pdf/grsvpol.f']
endif
if get_option('dssv')
    sources += ['Pdf/DSSV_gluon_update.f']
endif

# build
LeProHQpp = shared_module(
    'LeProHQpp',
    ['Python.cpp'] + sources,
    name_prefix: '',
    include_directories : [inc, configuration_inc],
    dependencies: [gsl_dep, gslcblas_dep, dvegas_dep, lhapdf_dep, boost_dep, rk_dep, py_dep, pye_dep, yaml_dep],
)