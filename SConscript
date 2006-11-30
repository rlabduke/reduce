Import("env_base", "env_etc")

reduce_scons_env = env_base.Copy(
  CXXFLAGS=env_etc.cxxflags_base,
  LIBS=env_etc.libm
)
if (env_etc.compiler == "unix_gcc"):
  reduce_scons_env.Append(CXXFLAGS=["-fno-strict-aliasing"])

Export("reduce_scons_env")

SConscript("toolclasses/SConscript")
SConscript("libpdb/SConscript")
SConscript("reduce_src/SConscript")
