# Ipopt Scilab toolbox

## Purpose 

The goal of this toolbox is to provide an interface to **IPOPT**,
an optimization solver based on an Interior Point algorithm.

This module is mainly based on compiled source code.

## Dependencies

In order to compile this module, you must compile its dependencies first (METIS, MUMPS, IPOPT). 

 * Linux and macOS

Go in the `thirdparty/build` directory and type
  
```
sh build.sh
```
If the build process is successfull, the libraries should be installed in `thirdparty/Linux/` or  
`thirdparty/Darwin/` depending on your platform

 * Windows
 
 For this operating system we use a prebuilt Ipopt library, which is available at
 
<https://github.com/coin-or/Ipopt/releases/download/releases%2F3.14.19/Ipopt-3.14.19-win64-msvs2022-md.zip>

Unzip this archive in the thirdparty folder and rename the obtained folder Ipopt-3.13.2-win64-msvs2019-md as "Windows".


## Build the toolbox

Once the dependencies are installed, execute `builder.sce` from within Scilab
and pray. If you want to package an Atoms module execute `makedist.sce`.

## Authors

2020-2025 - UTC - Stéphane MOTTELET
2010 - Consortium Scilab - Digiteo - Michael Baudin
2009-2010 - Consortium Scilab - Digiteo - Yann Collette
2008 - Yann Collette

## Licence

This toolbox is released under the GPL v3 licence.
