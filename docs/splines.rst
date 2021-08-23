.. _splines description:

=================
Model description
=================

1D Spline model
+++++++++++++++

The 1D spline model with coefficients :math:`S_{i, m}` with :math:`i` specifying the spline interval and :math:`m` (=0-3)
specifying the polynomial order is given by:

.. math::
    f_S(x)=\sum_{m=0}^3S_{i,m} \left(\frac{x-t_i}{\Delta t_i}\right)^m

    t_i\leq x \leq t_{i+1}, \Delta t_i=t_{i+1}-t_i

Each interval is represented by 4 coefficients :math:`S_{i,:}`. The positions :math:`t_i` specify the left borders
of the spline interval :math:`i` and :math:`\Delta t_i` is its size.

2D Spline model
+++++++++++++++

The 2D spline model with coefficients :math:`S_{i, j, m, n}` with :math:`(i,j)` specifying the spline interval in 2D and
:math:`(m, n)` (=0-3) specifying the polynomial orders in x and y is given by

.. math::
    f_S(x, y)=\sum_{m=0}^3\sum_{n=0}^3S_{i,j,m,n} \left(\frac{x-t_i}{\Delta t}\right)^m \left(\frac{y-u_j}{\Delta u}\right)^n

    t_i\leq x \leq t_{i+1}, \Delta t=t_{i+1}-t_i

    u_j\leq y \leq u_{j+1}, \Delta u=u_{j+1}-u_j

Each interval is represented by 16 coefficients :math:`S_{i,j,:,:}`. The positions :math:`(t_i, u_j)` specify the left borders
of the spline interval :math:`(i, j)` and :math:`(\Delta t_i, \Delta u_j)` is its size in x and y.

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