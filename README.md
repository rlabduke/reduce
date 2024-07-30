# reduce
Reduce - tool for adding and correcting hydrogens in PDB files

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

You can install reduce using conda package manager.

```shell
# Install latest version
conda install reduce -c bioconda
```

