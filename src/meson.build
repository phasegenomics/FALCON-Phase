project(
  'FALCON-Phase',
  'c',
  version : '0.1.0',
  default_options : [
    'buildtype=release',
    'c_std=gnu11',
   #'cpp_std=c++11',
    'b_ndebug=if-release',
    'warning_level=3',
  ],
  license : 'BSD',
  meson_version : '>= 0.41.1',
)

###########
# sources #
###########

sources = files([
    'bam_filt.c',
    'count_motif.c',
    'dna.c',
    'matrix.c',
    'run_phasing.c',
    'string_parser.c',
])

############
# CXXFLAGS #
############

warning_flags = []
compiler = meson.get_compiler('c')
foreach cflag: [
  #'-Wno-unused-result',
  '-Wno-unused-parameter',
  '-Wno-sign-compare',
  ]
    if compiler.has_argument(cflag)
      warning_flags += cflag
    endif
endforeach

if compiler.get_id() == 'clang'
  foreach cflag: [
    #'-Wno-unused-private-field',
    #'-Wno-unused-command-line-argument',
    #'-Wno-unused-variable',
    #'-Wno-unused-function',
    #'-Wno-char-subscripts',
    #'-Wno-pessimizing-move',
    #'-Wno-format-pedantic',
    #'-Wno-keyword-macro',
    ]
      if compiler.has_argument(cflag)
        warning_flags += cflag
      endif
  endforeach
endif

if compiler.get_id() == 'gcc'
  foreach cflag: [
    '-Wno-unused-but-set-variable',
    #'-Wno-unused-variable',
    #'-Wno-unused-function',
    #'-Wno-char-subscripts',
    #'-Wno-pedantic',
    #'-Wno-unused-but-set-variable',
    ]
      if compiler.has_argument(cflag)
        warning_flags += cflag
      endif
  endforeach
endif

compiler_flags = ['-march=native', '-D_GNU_SOURCE']
extra_flags = []

################
# Dependencies #
################

thread_dep = dependency('threads', required : true)
zlib_dep = dependency('zlib', required : true)
htslib_dep = dependency('htslib', required : true)
all_dependencies = [thread_dep, zlib_dep, htslib_dep]

###########
# Headers #
###########

#include_directories = [include_directories('.')]
include_directories = []

lib_install = (not meson.is_subproject()) or (get_option('default_library') == 'shared')

falcon_phase_lib = library(
  'falcon-phase',
  sources,
  soversion : 0,
  version : meson.project_version(),
  install : lib_install,
  dependencies : all_dependencies,
  include_directories : include_directories,
  c_args : [extra_flags, warning_flags, compiler_flags],
)
falcon_phase_bin = executable(
  'falcon-phase',
  sources + ['main.c'],
  install: true,
  dependencies : all_dependencies,
  include_directories : include_directories,
  link_with : [falcon_phase_lib],
  #override_options : ['c_std=c99'],
  c_args : [warning_flags, compiler_flags],
)

#################
# Maybe someday #
#################

# Google test.
#gtest_dep = dependency('gtest', main : true, required : false)
#if not gtest_dep.found()
#    gtest_proj = subproject('gtest')
#    gtest_inc = gtest_proj.get_variable('gtest_incdir')
#    gtest_lib = static_library('gtest', gtest_proj.get_variable('gtest_libsources'),
#                                gtest_proj.get_variable('gtest_mainsources'),
#                                include_directories : gtest_inc)
#
#    gtest_dep = declare_dependency(include_directories : gtest_inc,
#                               link_with : gtest_lib, dependencies: thread_dep)
#endif

# subdir('test_src')

# if gtest_dep.found()
#     test_localqv_bin = executable(
#         'test_falcon_phase',
#         sources + test_sources + ['test_src/main_test.cc'],
#         dependencies : all_dependencies,
#         include_directories : include_directories,
#         link_with : [falcon_phase_lib],
#         cpp_args : [warning_flags, compiler_flags])
# endif
