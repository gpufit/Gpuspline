.. _installation-and-testing:

========================
Installation and Testing
========================

Building from source code
+++++++++++++++++++++++++

This section describes how to build Gpuspline from source code. The source code has been build successfully and
been tested with the Microsoft Visual Studio compiler (2013 - 2019) on Windows 10 as well as with gcc 9.3 on Ubuntu 20.04.

Prerequisites
-------------

The following tools are required in order to build Gpuspline from source.

*Required*

* CMake_ 3.11 or later
* A C/C++ Compiler (gcc on Linux, Visual Studio on Windows)

*Optional*

* MATLAB_ if building the MATLAB bindings
* Python_ if building the Python bindings (and for producing Latex and HTML documentation output)
* PDF Latex installation (like Miktex on Windows or texlive-binaries on Linux) if converting the documentation from Latex to PDF

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
run cmake to configure and generate the compiler input files. 

Using the CMake Graphical User Interface
----------------------------------------

There is a graphical user interface available for CMake, which simplifies
the configuration and generation steps.  For further details, see
`Running CMake <https://cmake.org/runningcmake/>`_. The following steps outline 
how to use the basic features of the CMake GUI.

First, select the source code directory (the top level directory where the Gpuspline 
source code is located), and the build directory (where the binaries will be built).  
For this example, the source directory might be *C:\\Sources\\Gpuspline*, and the 
build directory might be *C:\\Sources\\Gpuspline-build-64*.

Next, click the "Configure" button, and select the desired compiler from the drop 
down list (e.g. Visual Studio 16 2019).  Under *Optional platform for Generator*, 
select the desired architecture (e.g. *x64* to compile 64-bit binaries).

Once configuration is complete, CMake will have automatically found the Matlab 
installation, and the installation directories will be listed in the *NAME* and 
*VALUE* columns.  If the Matlab installation was not found, the entries in the 
*VALUE* column can be manually edited.

Next, click on *Generate* to generate the Visual Studio solution files, which
will be used to build the Gpuspline package.

Running CMake from the command line
-----------------------------------

The following commands, executed from the command prompt, assume that the 
cmake executable (e.g. *C:\\Program Files\\CMake\\bin\\cmake.exe*) is automatically found
via the PATH environment variable (if not, the full path to cmake.exe must be
specified). This example also assumes that the source and build directories
have been set up as specified above.

.. code-block:: bash

    cd C:\Sources\Gpuspline-build-64
    cmake -G "Visual Studio 16 2019 Win64" C:\Sources\Gpuspline

Note that in this example the *-G* flag has been used to specify the
64-bit version of the Visual Studio 14 compiler. This flag should be changed
depending on the compiler used, and the desired architecture
(e.g. 32- or 64-bit). Further details of the CMake command line arguments
can be found `here <https://cmake.org/cmake/help/latest/manual/cmake.1.html>`__.


Common issues encountered during CMake configuration
----------------------------------------------------

**Python launcher**

Set Python_WORKING_DIRECTORY to a valid directory, it will be added to the
Python path.

**Matlab launcher**

Set Matlab_WORKING_DIRECTORY to a valid directory, it will be added to
the Matlab path.

Compiling Gpuspline on Windows
------------------------------

After configuring and generating the solution files using CMake, go to the
desired build directory and open Gpuspline.sln using Visual Studio. Select the
"Debug" or "Release" build options, as appropriate. Select the build target
"ALL_BUILD", and build this target. If the build process completes
without errors, the Gpuspline binary files will be created in the corresponding
"Debug" or "Release" folders in the build directory.

Compiling Gpuspline on Linux
----------------------------

The following commands can be executed to build Gpuspline on Linux.

.. code-block:: bash

    git clone https://github.com/gpufit/Gpuspline.git Gpuspline
    mkdir Gpuspline-build
    cd Gpuspline-build
    cmake -DCMAKE_BUILD_TYPE=RELEASE ../Gpuspline
    make

Run the test with

.. code-block:: bash

    ./splines_tests

To install the Python package

.. code-block:: bash

   cd pyGpuspline/dist
   pip install pyGpuspline-X.Y.Z-py2.py3-none-any.whl

Finally run the examples with (matplotlib needs a backend, for example pyqt5)

.. code-block:: bash

    pip install matplotlib
    pip install pyqt5
    python ../../../Gpuspline/examples/python/example_1d_interpolation.py
    python ../../../Gpuspline/examples/python/example_2d_resampling.py
   
Optional: Depending on the gcc and the Matlab versions, to run the Matlab package you may need to tell Matlab to use a
newer version of the C++ standard library

.. code-block:: bash

   export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6

Start Matlab.

.. code-block:: bash

   matlab

Then in Matlab add the matlab output directory to the path and execute some examples.

.. code-block:: bash

   addpath('XX/Gpuspline-build/matlab');
   cd('XX/Gpuspline/src/examples/matlab');
   example_1d_interpolation();
