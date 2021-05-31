project('LeProHQ',
    'cpp', 'fortran',
    version : '0.1.0',
    license : 'GPL3',
    default_options : 'cpp_std=c++11',
    meson_version : '>=0.55'
)

# mostly standard libraries
# https://github.com/GStreamer/gstreamer/blob/master/meson.build
cpp = meson.get_compiler('cpp')
gsl_dep = cpp.find_library('gsl')
gslcblas_dep = cpp.find_library('gslcblas')
boost_dep = dependency('boost', modules : ['system', 'filesystem', 'python3'])
yaml_dep = dependency('yaml-cpp')

# Python
py3_mod = import('python')
py3 = py3_mod.find_installation('python3')
py_dep = dependency('python3', version : '>=3.9')
pye_dep = dependency('python3-embed', version : '>=3.9')

# HEP libraries
lhapdf_dep = dependency('lhapdf', version : '>=6')
rk_dep = dependency('rk')
dvegas_dep = dependency('dvegas')

subdir('include')
subdir('src')
#subdir('tests')
