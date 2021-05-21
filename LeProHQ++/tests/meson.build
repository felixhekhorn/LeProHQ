
executable(
    'test_gsl',
    'gsl.cpp',
    include_directories : inc,
    dependencies: [gsl_dep, gslcblas_dep]
)

executable(
    'test_interpolation',
    ['interpolation.cpp', '../src/FullyInclusiveSoftResummed/Interpolation.cpp','../src/Common/EnvUtils.cpp','../src/Pdf/PdfWrapper.cpp'],
    include_directories : inc,
    dependencies: [gsl_dep, gslcblas_dep, lhapdf_dep, boost_dep]
)

executable(
    'test_me',
    ['me.cpp', '../src/FullyDiff/ME/A3.cpp'],
    include_directories : inc,
    dependencies: []
)

#executable(
#    'test_partonic',
#    ['partonic.cpp'] + sources,
#    include_directories : inc,
#    dependencies: []
#)