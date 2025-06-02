# Honeycomb2 <img src="resources/logo.png" alt="drawing" width="30"/>

This is the new version of the twist-3 evolution code [honeycomb](https://github.com/QCDatHT/honeycomb).
The new version sees a complete reworking of the discretization procedure and API. 
There is no backward compatibility. 

> The library is still work-in-progress. The API is somewhat stable, but major changes can still happen.

### Installation

**Prerequisite:** 
1. You should have installed a C++ compiler with support for C++20, including the `format` header.
2. You should have installed `cmake`

**Installation:**
The usual `cmake` procedure. From within `honeycomb2` directory
```shell
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make
make install
```
*Note:* By the default, if no prefix specification is given, `honeycomb2` will be installed in the `/usr/local/`  directory. If you want (or need) to use a different path, remember to export the `honeycomb2` `/lib` folder into the `LD_LIBRARY_PATH`.

After installation, you can run 
```shell
Honeycomb2-config --help
```
to get the list of available flags to be used in your project when using `honeycomb2`.

`honeycomb2` can be un-installed by running:
```shell
make clean
xargs rm < install_manifest.txt
```
from within the `build` directory.

**C/FORTRAN interface:**
`honeycomb2` supplies a limited API for C/FORTRAN programs. For C programs, one should 
```C
#include <honeycomb2/honeycomb2_c_api.h>
```
which expose the C-compatible API calls. 

In FORTRAN programs, one should remember to declare external the API calls, for instance
```FORTRAN
      external hc2_fi_set_up
      external hc2_fi_set_model
      external hc2_fi_evolve
      external hc2_fi_unload
      double precision hc2_fi_get_model
      external hc2_fi_get_model
```
where the trailing underscore is not present, since FORTRAN mingles function names automatically with the trailing underscore.

Both C and FORTRAN program must link with `libHoneycomb2.so`. The linker flags can be readily obtained running 
`Honeycomb2-config --ldflags`.
A minimalistic makefile for compiling a single FORTRAN file and link it with `honeycomb2` is:
```makefile
all: test

test: main.f
	gfortran $^ $(shell Honeycomb2-config --ldflags) -o $@

run: test
	LD_LIBRARY_PATH=$(shell Honeycomb2-config --libdir):$(LD_LIBRARY_PATH) ./test

```
### Examples 
Examples are under construction...
