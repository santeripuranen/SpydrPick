# SpydrPick

## Get the code
```
git clone --recursive https://github.com/santeripuranen/SpydrPick.git
```
See how to compile [here](README.md/#building-SpydrPick).


## About

Spydrpick is a tool for performing direct coupling analysis of aligned categorical datasets. It constructs a coupling graph based on thresholded node-pair (edge) scoring,
followed by pruning of edges that represent indirect couplings.


## Building SpydrPick

In order to compile the `SpydrPick` binary, go to the `build` directory (create one if necessary; in-source builds are strongly discouraged) and give these commands:

```
cmake ..
make SpydrPick
```

This should set up the CMake project, compile the binary and place it into the `bin` directory. If not, then take a look at [Compile-time dependencies](README.md/#compile-time-dependencies).

The `SpydrPick` binary will by default be statically linked, except for the standard C++ runtime library, which is unfeasible to link statically, and [TBB](https://www.threadingbuildingblocks.org/) that can only be linked dynamically. Installing SpydrPick to another location is as easy as copying the binary (given that the [TBB](https://www.threadingbuildingblocks.org/) runtime library is properly installed on your system).


### Compile-time dependencies

The SpydrPick code is written in C++ and wrapped into a [CMake](https://cmake.org/) project. It relies on several external libraries, most of which are fairly common in C++ software development. Your build environment must have the following requirements and compile-time dependencies satisfied:

* A C++14 compliant compiler (development was done using the [GNU C++ compiler](https://gcc.gnu.org/))
* [CMake](https://cmake.org/)
* [Boost](https://www.boost.org/)
* [Intel(R) Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org/)

You may need to set the [`CMAKE_MODULE_PATH`](https://cmake.org/cmake/help/latest/variable/CMAKE_MODULE_PATH.html) environment variable in order for CMake to find all relevant packages.


