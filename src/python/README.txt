Python binding for the [Gpuspline library](https://github.com/gpufit/Gpuspline) for the calculation of multidimensional cubic splines.

Installation

Currently the wheel file has to be installed locally.

If NumPy is not yet installed, install it using pip from the command line

pip install numpy

Then install pyGpufit from the local folder via:

pip install --no-index --find-links=LocalPathToWheelFile pyGpuspline

Examples

See /examples/python folder. Additionally require matplotlib and a matplotlib backend (tk or pyqt5).