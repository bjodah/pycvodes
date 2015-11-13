v0.2.3
======
- Better debugging (preserve back-trace from calling rhs() and jac())

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
