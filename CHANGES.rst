v0.6.0
======
- Bug-fix in get_integrator, dx_min and dx_max were ignored.
- Refactored to use AnyODE base class (share code with pyodeint & pygslodeiv2)

v0.5.0
======
- C++ wrapper API:
    - banded_padded_jac_cmaj -> banded_jac_cmaj
    - allow callbacks to indicate recoverable errors.

v0.4.4
======
- Better sdist

v0.4.3
======
- Better const correctness and other improvements in C++ wrapper

v0.4.2
======
- More robust setup.py

v0.4.1
======
- Added 'time_wall' output from integration.
- Added 'roots_output' to info dict of predefined

v0.4.0
======
- kwarg 'iterative' changed to 'iter_type' and 'linear_solver'
- sparse option dropped
- more flexible C++ interface
- pycvodes.get_include() useful for other projects linking against sundials (cvodes)

v0.3.0
======
- Better debugging (preserve back-trace from calling rhs() and jac())
- Changes to info dict: rename 'nrhs' -> 'nfev', 'njac' -> 'njev', added 'cpu_time', 'success'

v0.2.2
======
- Added support for root finding.
- Allow user to set maximum number of steps (previously only CVode's default of 500 was used).
- Improved derivative handling (for interpolation).
- Added option to make output from adaptive more sparse.

v0.2.1
======
- Added support for (first) derivative in output
- Min and max step now allowed to be set

v0.2.0
======
- New function signature: integrate_predefined and integrate_adaptive now
  also return an info dict containing ``nrhs`` and ``njac`` containing
  number of calls to each function made during last integration.
- Expose ``pycvodes.steppers`` tuple.
- check_callbable and check_indexing kwargs now defaults to False

v0.1.1
======
- Added lband, uband kwargs (compatible with scipy's signature)

v0.1
====
- Initial release
