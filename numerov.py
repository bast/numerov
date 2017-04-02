"""
Numerov: Compute vibrational levels, wavefunctions, and expectation values using the Numerov-Cooley algorithm.
"""

__all__ = [
    'solve_numerov',
    '__version__',
]

if __name__ == '__main__':
    # if run as a script or by 'python -m numerov'
    # we trigger the below "else" condition by the following import
    import numerov
    raise SystemExit(numerov.main())

# else we are imported
from _numerov.main import main
from _numerov.solve import solve_numerov

from _numerov import __version__
