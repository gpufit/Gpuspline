# Creating a binary package/release

The binary package bundles different builds outputs into a single distributable binary package containing the Gpuspline library,
the Matlab and Python bindings and the documentation.

Follow this step by step recipe to create a Windows binary package.

## Set/Update the version number

Unfortunately the version has to be updated in various places.

- CMakeLists.txt (project( Gpuspline VERSION 1.0.0 ))
- docs/conf.py (release = u'1.0.0')
- src/matlab/spline_version.m 
- src/python/pygpuspline/version.py
- in the call to the packaging script (create_package.bat %1 1.0.0 %3)
- package/sdk_readme.txt

Push to Github afterwards (you can add a Git tag for the new version with git tag -a vX.Y.Z -m "Version X.Y.Z release" followed by git push --tags).

## Convert Documentation from restructured text to html/latex

Use sphinx and docs/make.bat with parameters html and latex.

## Use CMAKE to generate the project

- Build directory for Visual Studio 2019 Win64 is BUILD_BASE_PATH/VC16x64
- Matlab, Python, Latex (e.g. Miktex) must be available

See also [Build from sources](https://gpuspline.readthedocs.io/en/latest/installation.html#building-from-source-code) for instructions.

## Build for Win64

Build with single precision. Everything should succeed.

- Configuration RelWithDebInfo is used for all builds!
- With Visual Studio 2019 Win64 build target PYTHON_WHEEL, MATLAB_SPLINE_PACKAGE and DOCUMENTATION_PDFLATEX (SOURCE_BASE_PATH\docs\_build\latex\Gpuspline.pdf will be created from Gpuspline.tex at the same location)

## Run the examples for the Bindings

In Matlab and Python.

## Call the assemble script

create_package.bat %1 %2 %3

with 

- %1 is the BUILD_BASE_PATH (the path containing the various (see below) CMake generated Visual Studio projects)

- %2 is the VERSION (e.g. 1.0.0)

- %3 is the SOURCE_BASE_PATH (the path containing the sources)

The output is a folder (BUILD_BASE_PATH/Gpuspline-VERSION) which is also zipped if 7-Zip is available.

## Retrieve the hash for the current commit in GIT

git rev-parse --verify --short HEAD

and add as suffix to build file name (this specifies the source commit state)