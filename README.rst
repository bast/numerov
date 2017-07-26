

Numerov
=======

.. image:: https://travis-ci.org/bast/numerov.svg?branch=master
   :target: https://travis-ci.org/bast/numerov/builds

.. image:: https://coveralls.io/repos/github/bast/numerov/badge.svg?branch=master
   :target: https://coveralls.io/github/bast/numerov?branch=master

.. image:: https://img.shields.io/badge/license-%20MPL--v2.0-blue.svg
   :target: https://github.com/bast/numerov/blob/master/LICENSE

Compute vibrational levels, wavefunctions, and expectation values using the
Numerov-Cooley algorithm.


Copyright and license
---------------------

Copyright 2017 Radovan Bast.
Use of this source code is governed by a the Mozilla Public License v2.0 that
can be found in the
`LICENSE file <https://github.com/bast/numerov/blob/master/LICENSE>`_.


Citation
--------

If you use this tool in a program or publication, please acknowledge its
author(s) by adding the following reference (please replace vX.Y.Z
by the appropriate version):

- Radovan Bast, Numerov vX.Y.Z, 2017, GitHub repository, https://github.com/bast/numerov.

As soon as the API stabilizes a bit, I will tag a release and attach a DOI to it.


Installation
------------

.. code-block:: bash

   pip install numerov


Example
-------

Please have a look `here <https://github.com/bast/numerov/tree/master/pnc-example>`_.


Background
----------

This script will calculate the vibrational levels (and wavefunctions)
corresponding to a normal mode numerically using the Numerov-Cooley algorithm.

The script will increase the energy and count the nodes of the wave function.
If the number of nodes changes and stepsize is below ``energy_precision_hartree``, it will
accept the solution, integrate the property along q and move on to the next
solution until ``num_solutions`` is reached. It will also calculate the transition
frequencies 0 -> n, this is useful to check against the harmonic frequencies.


Advice
------

- ``energy_precision_hartree`` is often more important than number of grid points.
- Be careful with the displacement range.
- If the script enters an endless loop probably the reduced mass or the displacement range is wrong.
- Practice first with the harmonic oscillator.
- It is a good idea to play with parameters to check convergence and numerical
  stability.
- Potential and property are approximated by polynomials that contain
  coefficients FROM ZEROTH to nth order (that's what polyfit gives) and you
  might not want that (for instance you might insist that the gradient of the
  potential is zero at equilibrium) in this case you can provide your own
  expansion coefficients or program an alternative interpolation scheme.
