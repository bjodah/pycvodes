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
