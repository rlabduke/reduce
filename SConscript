import libtbx
Import("env_base", "env_etc")

reduce_scons_env = env_base.Clone(
  LIBS=env_etc.libm)
if (libtbx.manual_date_stamp < 20090819):
  # XXX backward compatibility 2009-08-19
  env.Replace(CCFLAGS=env_etc.ccflags_base)
  env.Replace(CXXFLAGS=env_etc.cxxflags_base)

Export("reduce_scons_env")

SConscript("toolclasses/SConscript")
SConscript("libpdb/SConscript")
SConscript("reduce_src/SConscript")
