.. _installation-and-testing:

========================
Installation and Testing
========================

Building from source code
+++++++++++++++++++++++++

This section describes how to build Gpuspline from source code. Note that as of
the initial release of Gpuspline, the source code has been tested only with the
Microsoft Visual Studio compiler.

Prerequisites
-------------

The following tools are required in order to build Gpuspline from source.

*Required*

* CMake_ 3.11 or later
* A C/C++ Compiler

  * Linux: ...
  * Windows: Visual Studio 2013, 2015 or 2017

*Optional*

* MATLAB_ if building the MATLAB bindings (minimum version Matlab ...)
* Python_ if building the Python bindings (Python version 2.x or 3.x)
* PDF Latex installation (like Miktex) if converting the documentation from Latex to PDF

Source code availability
------------------------

The source code is available in an open repository hosted at Github, at the
following URL.

.. code-block:: bash

    https://github.com/gpufit/Gpuspline.git

To obtain the code, Git may be used to clone the repository.

Compiler configuration via CMake
--------------------------------

CMake is an open-source tool designed to build, test, and package software.
It is used to control the software compilation process using compiler
independent configuration files, and generate native makefiles and workspaces
that can be used in the compiler environment. In this section we provide a
simple example of how to use CMake in order to generate the input files for the
compiler (e.g. the Visual Studio solution file), which can then be used to
compile Gpuspline.

First, identify the directory which contains the Gpuspline source code
(for example, on a Windows computer the Gpuspline source code may be stored in
*C:\\Sources\\Gpuspline*). Next, create a build directory outside the
source code source directory (e.g. *C:\\Sources\\Gpuspline-build-64*). Finally,
run cmake to configure and generate the compiler input files. The following
commands, executed from the command prompt, assume that the cmake executable
(e.g. *C:\\Program Files\\CMake\\bin\\cmake.exe*) is automatically found
via the PATH environment variable (if not, the full path to cmake.exe must be
specified). This example also assumes that the source and build directories
have been set up as specified above.

.. code-block:: bash

    cd C:\Sources\Gpuspline-build-64
    cmake -G "Visual Studio 14 2015 Win64" C:\Sources\Gpuspline

Note that in this example the *-G* flag has been used to specify the
64-bit version of the Visual Studio 14 compiler. This flag should be changed
depending on the compiler used, and the desired architecture
(e.g. 32- or 64-bit). Further details of the CMake command line arguments
can be found `here <https://cmake.org/cmake/help/latest/manual/cmake.1.html>`__.

There is also a graphical user interface available for CMake, which simplifies
the configuration and generation steps. For further details, see
`Running CMake <https://cmake.org/runningcmake/>`_.

Common issues encountered during CMake configuration
----------------------------------------------------

**Python launcher**

Set Python_WORKING_DIRECTORY to a valid directory, it will be added to the
Python path.

**Matlab launcher**

Set Matlab_WORKING_DIRECTORY to a valid directory, it will be added to
the Matlab path.

Compiling Gpuspline on Windows
------------------------------------

After configuring and generating the solution files using CMake, go to the
desired build directory and open Gpuspline.sln using Visual Studio. Select the
"Debug" or "Release" build options, as appropriate. Select the build target
"ALL_BUILD", and build this target. If the build process completes
without errors, the Gpuspline binary files will be created in the corresponding
"Debug" or "Release" folders in the build directory.

Compiling Gpuspline on Linux
----------------------------------
 ...

MacOS
-----

Gpuspline has not yet been officially tested on a computer running MacOS.
However, satisfying the Prerequisites_ and using CMake, we estimate that the
library should build in principle and one should also be able to run the
examples on MacOS.
