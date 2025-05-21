# reduce
Reduce - tool for adding and correcting hydrogens in PDB files

## Deprecation notice

This version of reduce is being replaced by reduce2, which is integrated into the CCTBX
project as described at https://github.com/cctbx/cctbx_project/tree/master/mmtbx/reduce
with additional development being done there.  The reduce2 program handle native CIF file
input and is being maintained with ongoing file-format and other changes.  Reduce2 is
slower than reduce and there are a few corner-case files where it does not work, but it
produces more accurate results in most cases.

## Building

### Makefile

The original Makefile-based build can be used on Unix-like operating systems.  To use it, you run `make` in
the root directory.  It will build libraries in the **libtbx**, **toolclasses**, and **reduce_src** directories and
place the **reduce** executable in the **reduce_src** subdirectory.

Once the code has been installed, run `sudo make install` to install it in /usr/local/bin and the
required files in /usr/local.

### CMake

A CMakeLists.txt file has been created to support cross-platform building use CMake.  This supports
building on Windows in addition to Unix-like platforms.

This approach supports out-of-source builds, where the build is done in a directory other than
the source directory, so that build artifacts do not clutter the source tree (and in particular,
to avoid clobbering the original Makefile files).

If the source has been checked out in **~/src/reduce** then the following is one approach to doing
a CMake-based build:

    mkdir -p ~/build/reduce
    cd ~/build/reduce
    cmake ~/src/reduce
    make
    sudo make install

This will install the executable by default in /usr/local/bin and the required files in /usr/local.
As with other CMake projects, the you can set the **CMAKE_INSTALL_PREFIX** variable and it will put the
executable there under the bin/ directory and the HET dictionary in the prefix directory,
telling the program where to find it.

### Running

Once installed, if /usr/local/bin is on your PATH, reduce can be run using `reduce`.
Use `reduce -help` for more information about the arguments and behavior of the program.

## Installation using package manager

Ueno, M. (github ID eunos-1128) has graciously made the original (now deprecated) version of
reduce available through standard installation channels as described below.

### Conda

You can install reduce on Linux and Mac using the conda package manager.  On Windows, you can
use the Windows Subsystem for Linux and install reduce using conda in the WSL environment.

```shell
conda install reduce -c bioconda
```

### Homebrew/Linuxbrew
You can also install reduce using homebrew (or linuxbrew).

```shell
brew tap brewsci/bio
brew install reduce
```
