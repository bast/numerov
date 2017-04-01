[![License](https://img.shields.io/badge/license-%20MPL--v2.0-blue.svg)](../master/LICENSE)


# Numerov

Compute vibrational levels (and wavefunctions) using the Numerov-Cooley algorithm.

This script will calculate the vibrational levels (and wavefunctions)
corresponding to a normal mode numerically using the Numerov-Cooley algorithm.

The script will increase the energy and count the nodes of the wave function.
If the number of nodes changes and stepsize is below energy_precision, it will
accept the solution, integrate the property along q and move on to the next
solution until nr_solutions is reached.  It will also calculate the transition
frequencies 0 -> n this is useful to check against the harmonic frequencies.

Advice:

- Energy_precision is often more important than number of grid points.
- Be careful with the range (q_min and q_max).
- If the script enters an endless loop probably the reduced mass is wrong or
  the range (q_min and q_max) is too large.
- Practice first with the harmonic oscillator.
- It is a good idea to play with parameters to check convergence and numerical
  stability.
- Potential and property are approximated by polynomials that contain
  coefficients FROM ZEROTH to nth order (that's what polyfit gives) and you
  might not want that (for instance you might insist that the gradient of the
  potential is zero at equilibrium) in this case you can provide your own
  expansion coefficients or program an alternative interpolation scheme.
