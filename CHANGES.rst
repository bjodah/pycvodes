v0.14.0
=======
- Allow use of jtimes_setup from CVodes (currently only from C++).

v0.13.1
=======
- Fix in bundled AnyODE

v0.13.0
=======
- Reorganized pycvodes.config

v0.12.2
=======
- Python compatibility fix

v0.12.1
=======
- Update setup.py

v0.12.0
=======
- Require sundials 5.1.0

v0.11.15
========
- stab_lim_det

v0.11.13
========
- Expose set_max_steps_between_jac

v0.11.12
========
- Fix parallel adaptive wrt ew_ele array (triple pointer -> double pointer).

v0.11.11
========
- Add compile time option "PYCVODES_CLIP_TO_CONSTRAINTS" (currently defaults to 0).
  To handle polynomial interpolatino from sundials giving negative values even when
  using positive constraints.

v0.11.10
========
- Support for sundials>=4.0.0

v0.11.9
=======
- Fix signature in header

v0.11.8
=======
- Fix syntax error in header file

v0.11.7
=======
- Add support for constraints when using sundials>=3.2.0

v0.11.6
=======
- Bug fix in adaptive (gh-92 by @spcornelius)

v0.11.5
=======
- Add REAL_TYPE, INDEX_TYPE variables to _config.py

v0.11.4
=======
- Bug fix (const qualifier, cf. https://github.com/bjodah/pycvodes/pull/84)

v0.11.3
=======
- New python signature: t is now a NumPy scalar

v0.11.2
=======
- New AnyODE version

v0.11.1
=======
- Macro USE_LAPACK renamed (with inverted logic) to PYCVODES_NO_LAPACK
- New Macro: PYCVODES_NO_KLU (for compiling with LAPACK but without KLU)

v0.11.0
=======
- Enhanced _config.py
- Support for sparse jacobian

v0.10.12
========
- Support for constraints in sundials 3.2.0

v0.10.11
========
- Fix to setup.py

v0.10.10
========
- Determine compilation options by attempt in sdist

v0.10.9
=======
- More robust ``pycvodes.config``.

v0.10.8
=======
- Bug fix in adaptive (gh-93 by @spcornelius)

v0.10.7
=======
- LAPACK now optional
- Builds on Windows

v0.10.6
=======
- Optionally return error weights and estimated local errors (kwarg: ``ew_ele=True``)

v0.10.5
=======
- Bump to AnyODE 14
- Only require C++14

v0.10.4
=======
- Build related changes.

v0.10.3
=======
- More robust deducation of sundials version.

v0.10.2
=======
- New AnyODE version (13)

v0.10.1
=======
- Better compile-time inspection of sundials version
- Work aroud for sundials/lapack inconsistency wrt. dgbtrf

v0.10.0
=======
- Exposed quadrature integration (quads and get_nquads)
- Bumped AnyODE version (12)
- More timing data (time spent in rhs, jac & preconditioners)

v0.9.2
======
- Bump AnyODE version

v0.9.1
======
- variable tidx exposed in simple_adaptive

v0.9.0
======
- adaptive integration now reallocs its own space (allows direct transfer of ownership to e.g. numpy arrays)

v0.8.4
======
- Setting the environment variable ANYODE_VERBOSITY to 0 now silences errors & warnings.

v0.8.3
======
- Add jtimes=False default option in simple_{adaptive,predefined}

v0.8.2
======
- Added ``record_steps`` option.

v0.8.1
======
- Explicit use of std::make_unique from the C++14 standard.

v0.8.0
======
- Use new (templated) AnyODE.
- Fix back-stepping logic in adaptive.

v0.7.6
======
- return nreached in parallel predefined

v0.7.5
======
- Return atol & rtol in info dict
- Fix 'success' in info dict when return_on_error & return_on_root are both true.

v0.7.4
======
- Add return_on_error to cvodes_anyode_parallel
- Use environment variable ANYODE_NUM_THREADS

v0.7.3
======
- support for record_rhs_xvals/record_jac_xvals/record_order/record_fpe

v0.7.2
======
- Address VisibleDeprecationWarning from numpy ndarray.reshape

v0.7.1
======
- get_dx_max_cb (callback to calculate dx_max)

v0.7.0
======
- dx0cb
- atol may now be vector even from Python

v0.6.1
======
- New kwarg for autonomous systems: autorestart=<int>, helps when h/t ~= machine epsilon
- New kwarg for ``adaptive``: return_on_error, useful to take a fixed number of steps.
- New non-public module: _config (stores choice of lapack for now)
- adaptive in cvodes_cxx now return starting point when x0 >= xend (was CV_ILL_INPUT)

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
