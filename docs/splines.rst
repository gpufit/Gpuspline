.. _splines description:

=================
Model description
=================

From a given set of d-dimensional (d=1-3) data points, a cubic spline representation is calculated with 4^d
coefficients per d-dimensional data interval. The spline model function will take exactly the value of the data points
and the end condition is that 1st and 2nd derivative of the spline model function is zero at the end points (and outside
the spline will be extrapolated by constant values, i.e. the nearest data values).
As reference the `Wolfram Mathworld page on Cubic splines <https://mathworld.wolfram.com/CubicSpline.html>`_ as well
as the `Python implementation of cubic splines <https://github.com/ZhuangLab/storm-analysis/tree/master/storm_analysis/spliner>`_ by the ZhuangLab were used.

The cubic spline representation uses the logical spacing of the data (i.e. the interval lengths are all one) with indexing
starting at 0. Please take care to translate the spline interval coordinates to your real x,y,z coordinates yourself.
In the following the spline model functions for d=1-3 dimensions is given explicitly. The spline model is fully
characterized by the set of spline coefficients :math:`S_{..}`.

1D Spline model
+++++++++++++++

The 1D spline model with coefficients :math:`S_{i, m}` with :math:`i` specifying the spline interval and :math:`m` (=0-3)
specifying the polynomial order is given by:

.. math::
    f_S(x)=\sum_{m=0}^3S_{i,m} \left(\frac{x-t_i}{\Delta t_i}\right)^m

    t_i\leq x \leq t_{i+1}, \Delta t_i=t_{i+1}-t_i

Each interval is represented by 4 coefficients :math:`S_{i,:}`. The positions :math:`t_i` specify the left borders
of the spline interval :math:`i` and :math:`\Delta t_i` is its size. Here :math:`t_i` goes from 0 to N-1 for N data points
and :math:`\Delta t_i` equals 1.

2D Spline model
+++++++++++++++

The 2D spline model with coefficients :math:`S_{i, j, m, n}` with :math:`(i,j)` specifying the spline interval in 2D and
:math:`(m, n)` (=0-3) specifying the polynomial orders in x and y is given by

.. math::
    f_S(x, y)=\sum_{m=0}^3\sum_{n=0}^3S_{i,j,m,n} \left(\frac{x-t_i}{\Delta t}\right)^m \left(\frac{y-u_j}{\Delta u}\right)^n

    t_i\leq x \leq t_{i+1}, \Delta t=t_{i+1}-t_i

    u_j\leq y \leq u_{j+1}, \Delta u=u_{j+1}-u_j

Each interval is represented by 16 coefficients :math:`S_{i,j,:,:}`. The positions :math:`(t_i, u_j)` specify the left borders
of the spline interval :math:`(i, j)` and :math:`(\Delta t_i, \Delta u_j)` is its size in x and y. Here :math:`(t_i, u_j)` goes from (0,0) to (N-1, M-1) for NxM data points
and :math:`(\Delta t_i, \Delta u_j)` equals (1,1).

3D Spline model
+++++++++++++++

The 3D spline model with coefficients :math:`S_{i, j, k, m, n, o}` with :math:`(i, j, k)` specifying the spline interval
in 3D and :math:`(m, n, o)` (=0-3) specifying the polynomial orders in x, y and z is given by

.. math::
    f_S(x, y, z)=\sum_{m=0}^3\sum_{n=0}^3\sum_{o=0}^3S_{i,j,k,m,n, o} \left(\frac{x-t_i}{\Delta t}\right)^m \left(\frac{y-u_j}{\Delta u}\right)^n \left(\frac{z-v_k}{\Delta v}\right)^o

    t_i\leq x \leq t_{i+1}, \Delta t=t_{i+1}-t_i

    u_j\leq y \leq u_{j+1}, \Delta u=u_{j+1}-u_j

    v_k\leq z \leq v_{k+1}, \Delta v=v_{k+1}-v_k

Each interval is represented by 64 coefficients :math:`S_{i,j,k,:,:,:}`. The positions :math:`(t_i, u_j, v_k)` specify the left borders
of the spline interval :math:`(i, j, k)` and :math:`(\Delta t_i, \Delta u_j, \Delta v_k)` is its size in x, y and z.
Here :math:`(t_i, u_j, v_k)` goes from (0,0,0) to (N-1, M-1, K-1) for NxMxK data points
and :math:`(\Delta t_i, \Delta u_j, \Delta v_k)` equals (1,1,1).