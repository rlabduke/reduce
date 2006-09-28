import os

reduce_scons_env = Environment(
  ENV=os.environ,
  tools=["cc", "g++", "gnulink", "ar"])

Export("reduce_scons_env")

SConscript("toolclasses/SConscript")
SConscript("libpdb/SConscript")
SConscript("reduce_src/SConscript")
