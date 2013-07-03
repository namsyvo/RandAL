def build(bld):
  bld(features     ='cxx cshlib',
      source       = 'wat_array.cpp bit_array.cpp fmIndex.cpp',
      name         = 'fmIndex',
      target       = 'fmIndex',
      includes     = '.')
  bld(features     ='cxx cstaticlib',
      source       = 'wat_array.cpp bit_array.cpp fmIndex.cpp',
      name         = 'wat_array',
      target       = 'wat_array',
      includes     = '.')
  bld(features     ='cxx cprogram',
      source       = 'wat_array.cpp bit_array.cpp fmIndex.cpp construct.cpp',
      name         = 'fmconstruct',
      target       = 'fmconstruct',
      includes     = '.')
  bld(features     ='cxx cprogram',
      source       = 'wat_array.cpp bit_array.cpp fmIndex.cpp search.cpp',
      name         = 'fmsearch',
      target       = 'fmsearch',
      includes     = '.')
  
  bld.install_files('${PREFIX}/include/fmIndex++', bld.path.ant_glob('*.hpp'))
  bld.install_files('${PREFIX}/usr/bin/', bld.path.ant_glob('construct'))
  bld.install_files('${PREFIX}/usr/bin/', bld.path.ant_glob('search'))
