# SUNDIALS Changelog

## Changes to SUNDIALS in release 7.4.0

### New Features and Enhancements

`ARKodeSetCFLFraction` now allows `cfl_frac` to be greater than or equal to one.

Added an option to enable compensated summation of the time accumulator for all
of ARKODE. This was previously only an option for the SPRKStep module. The new
function to call to enable this is `ARKodeSetUseCompensatedSums`.

### Bug Fixes

Fixed segfaults in `CVodeAdjInit` and `IDAAdjInit` when called after adjoint
memory has been freed.

Fixed a CMake bug that would cause the Caliper compile test to fail at configure
time.

Fixed a bug in the CVODE/CVODES `CVodeSetEtaFixedStepBounds` function which
disallowed setting `eta_min_fx` or `eta_max_fx` to 1.

`SUNAdjointStepper_PrintAllStats` was reporting the wrong quantity for the
number of "recompute passes" and has been fixed.

### Deprecation Notices

The `SPRKStepSetUseCompensatedSums` function has been deprecated. Use the
`ARKodeSetUseCompensatedSums` function instead.

## Changes to SUNDIALS in release 7.3.0

### Major Features

A new discrete adjoint capability for explicit Runge-Kutta methods has been
added to the ARKODE ERKStep and ARKStep stepper modules. This is based on a
new set of shared classes, `SUNAdjointStepper` and `SUNAdjointCheckpointScheme`.
A new example demonstrating this capability can be found in
`examples/arkode/C_serial/ark_lotka_volterra_ASA.c`.

### New Features and Enhancements

#### ARKODE

The following changes have been made to the default ERK, DIRK, and ARK methods
in ARKODE to utilize more efficient methods:

| Type               | Old Default                                                      | New Default                                                        |
| ------------------ | ---------------------------------------------------------------- | ------------------------------------------------------------------ |
| 2nd Order Explicit | `ARKODE_HEUN_EULER_2_1_2`                                        | `ARKODE_RALSTON_3_1_2`                                             |
| 4th Order Explicit | `ARKODE_ZONNEVELD_5_3_4`                                         | `ARKODE_SOFRONIOU_SPALETTA_5_3_4`                                  |
| 5th Order Explicit | `ARKODE_CASH_KARP_6_4_5`                                         | `ARKODE_TSITOURAS_7_4_5`                                           |
| 6th Order Explicit | `ARKODE_VERNER_8_5_6`                                            | `ARKODE_VERNER_9_5_6`                                              |
| 8th Order Explicit | `ARKODE_FEHLBERG_13_7_8`                                         | `ARKODE_VERNER_13_7_8`                                             |
| 2nd Order Implicit | `ARKODE_SDIRK_2_1_2`                                             | `ARKODE_ARK2_DIRK_3_1_2`                                           |
| 3rd Order Implicit | `ARKODE_ARK324L2SA_DIRK_4_2_3`                                   | `ARKODE_ESDIRK325L2SA_5_2_3`                                       |
| 4th Order Implicit | `ARKODE_SDIRK_5_3_4`                                             | `ARKODE_ESDIRK436L2SA_6_3_4`                                       |
| 5th Order Implicit | `ARKODE_ARK548L2SA_DIRK_8_4_5`                                   | `ARKODE_ESDIRK547L2SA2_7_4_5`                                      |
| 4th Order ARK      | `ARKODE_ARK436L2SA_ERK_6_3_4` and `ARKODE_ARK436L2SA_DIRK_6_3_4` | `ARKODE_ARK437L2SA_ERK_7_3_4` and `ARKODE_ARK437L2SA_DIRK_7_3_4`   |
| 5th Order ARK      | `ARKODE_ARK548L2SA_ERK_8_4_5` and `ARKODE_ARK548L2SA_DIRK_8_4_5` | `ARKODE_ARK548L2SAb_ERK_8_4_5` and `ARKODE_ARK548L2SAb_DIRK_8_4_5` |


The old default methods can be loaded using the functions `ERKStepSetTableName`
or `ERKStepSetTableNum` with ERKStep and `ARKStepSetTableName` or
`ARKStepSetTableNum` with ARKStep and passing the desired method name string or
constant, respectively. For example, the following call can be used to load the
old default fourth order method with ERKStep:
```
/* Load the old 4th order ERK method using the table name */
ierr = ERKStepSetTableName(arkode_mem, "ARKODE_ZONNEVELD_5_3_4");
```
Similarly with ARKStep, the following calls can be used for ERK, DIRK, or ARK
methods, respectively:
```
/* Load the old 4th order ERK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_DIRK_NONE",
                           "ARKODE_ZONNEVELD_5_3_4");

/* Load the old 4th order DIRK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_SDIRK_5_3_4",
                           "ARKODE_ERK_NONE");

/* Load the old 4th order ARK method by name */
ierr = ARKStepSetTableName(arkode_mem, "ARKODE_ARK436L2SA_DIRK_6_3_4",
                           "ARKODE_ARK436L2SA_ERK_6_3_4");
```

Additionally, the following changes have been made to the default time step
adaptivity parameters in ARKODE:

| Parameter             | Old Default           | New Default |
|-----------------------|-----------------------|-------------|
| Controller            | PID (PI for ERKStep)  | I           |
| Safety Factor         | 0.96                  | 0.9         |
| Bias                  | 1.5 (1.2 for ERKStep) | 1.0         |
| Fixed Step Bounds     | [1.0, 1.5]            | [1.0, 1.0]  |
| Adaptivity Adjustment | -1                    | 0           |

The following calls can be used to restore the old defaults for ERKStep:
```
SUNAdaptController controller = SUNAdaptController_Soderlind(ctx);
SUNAdaptController_SetParams_PI(controller, 0.8, -0.31);
ARKodeSetAdaptController(arkode_mem, controller);
SUNAdaptController_SetErrorBias(controller, 1.2);
ARKodeSetSafetyFactor(arkode_mem, 0.96);
ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
```
The following calls can be used to restore the old defaults for other ARKODE
integrators:
```
SUNAdaptController controller = SUNAdaptController_PID(ctx);
ARKodeSetAdaptController(arkode_mem, controller);
SUNAdaptController_SetErrorBias(controller, 1.5);
ARKodeSetSafetyFactor(arkode_mem, 0.96);
ARKodeSetFixedStepBounds(arkode_mem, 1, 1.5);
ARKodeSetAdaptivityAdjustment(arkode_mem, -1);
```
In both cases above, destroy the controller at the end of the run with
`SUNAdaptController_Destroy(controller);`.

The Soderlind time step adaptivity controller now starts with an I controller
until there is sufficient history of past time steps and errors.

Added `ARKodeSetAdaptControllerByName` to set a time step adaptivity controller
with a string. There are also four new controllers: `SUNAdaptController_H0211`,
`SUNAdaptController_H0321`, `SUNAdaptController_H211`, and
`SUNAdaptController_H312`.

Added the `ARKODE_RALSTON_3_1_2` and `ARKODE_TSITOURAS_7_4_5` explicit
Runge-Kutta Butcher tables.

Improved the precision of the coefficients for `ARKODE_ARK324L2SA_ERK_4_2_3`,
`ARKODE_VERNER_9_5_6`, `ARKODE_VERNER_10_6_7`, `ARKODE_VERNER_13_7_8`,
`ARKODE_ARK324L2SA_DIRK_4_2_3`, and `ARKODE_ESDIRK324L2SA_4_2_3`.

#### CVODE / CVODES

Added support for resizing CVODE and CVODES when solving initial value problems
where the number of equations and unknowns changes over time. Resizing requires
a user supplied history of solution and right-hand side values at the new
problem size, see `CVodeResizeHistory` for more information.

#### KINSOL

Added support in KINSOL for setting user-supplied functions to compute the
damping factor and, when using Anderson acceleration, the depth in fixed-point
or Picard iterations. See `KINSetDampingFn` and `KINSetDepthFn`, respectively,
for more information.

#### SUNDIALS Types

A new type, `suncountertype`, was added for the integer type used for counter
variables. It is currently an alias for `long int`.

### Bug Fixes

#### ARKODE

Fixed bug in `ARKodeResize` which caused it return an error for MRI methods.

Removed error floors from the `SUNAdaptController` implementations which could
unnecessarily limit the time size growth, particularly after the first step.

Fixed bug in `ARKodeSetFixedStep` where it could return `ARK_SUCCESS` despite
an error occurring.

Fixed bug in the ARKODE SPRKStep `SPRKStepReInit` function and `ARKodeReset` function
with SPRKStep that could cause a segmentation fault when compensated summation is not
used.

#### KINSOL

Fixed a bug in KINSOL where an incorrect damping parameter is applied on the
initial iteration with Anderson acceleration unless `KINSetDamping` and
`KINSetDampingAA` are both called with the same value when enabling damping.

Fixed a bug in KINSOL where errors that occurred when computing Anderson
acceleration were not captured.

Added missing return values to `KINGetReturnFlagName`.

#### CMake

Fixed the behavior of `SUNDIALS_ENABLE_ERROR_CHECKS` so additional runtime error
checks are disabled by default with all release build types. Previously,
`MinSizeRel` builds enabled additional error checking by default.

### Deprecation Notices

All work space functions, e.g., `CVodeGetWorkSpace` and `ARKodeGetLinWorkSpace`,
have been deprecated and will be removed in version 8.0.0.

## Changes to SUNDIALS in release 7.2.1

### New Features and Enhancements

Unit tests were separated from examples. To that end, the following directories
were moved out of the `examples/` directory to the `test/unit_tests` directory:
`nvector`, `sunmatrix`, `sunlinsol`, and `sunnonlinsol`.

### Bug Fixes

Fixed a bug in ARKStep where an extra right-hand side evaluation would occur
each time step when enabling the ``ARKodeSetAutonomous`` option and using an
IMEX method where the DIRK table has an implicit first stage and is not stiffly
accurate.

## Changes to SUNDIALS in release 7.2.0

### Major Features

Added a time-stepping module to ARKODE for low storage Runge-Kutta methods,
LSRKStep.  This currently supports five explicit low-storage methods: the
second-order Runge-Kutta-Chebyshev and Runge-Kutta-Legendre methods, and the
second- through fourth-order optimal strong stability preserving Runge-Kutta
methods. All methods include embeddings for temporal adaptivity.

Added an operator splitting module, SplittingStep, and forcing method module,
ForcingStep, to ARKODE. These modules support a broad range of operator-split
time integration methods for multiphysics applications.

Added support for multirate time step adaptivity controllers, based on the
recently introduced `SUNAdaptController` base class, to ARKODE's MRIStep module.
As a part of this, we added embeddings for existing MRI-GARK methods, as well as
support for embedded MERK and IMEX-MRI-SR methods. Added new default MRI methods
for temporally adaptive versus fixed-step runs.

### New Features and Enhancements

#### Logging

The information level logging output in ARKODE, CVODE(S), and IDA(S) has been
updated to be more uniform across the packages and a new `tools` directory has
been added with a Python module, `suntools`, containing utilities for parsing
logging output. The Python utilities for parsing CSV output have been relocated
from the `scripts` directory to the Python module.

#### SUNStepper

Added the `SUNStepper` base class to represent a generic solution procedure for
IVPs. This is used by the SplittingStep and ForcingStep modules of ARKODE. A
SUNStepper can be created from an ARKODE memory block with the new function
`ARKodeCreateSUNStepper`. To enable interoperability with `MRIStepInnerStepper`,
the function `MRIStepInnerStepper_CreateFromSUNStepper` was added.

#### ARKODE

Added functionality to ARKODE to accumulate a temporal error estimate over
multiple time steps. See the routines `ARKodeSetAccumulatedErrorType`,
`ARKodeResetAccumulatedError`, and `ARKodeGetAccumulatedError` for details.

Added the `ARKodeSetStepDirection` and `ARKodeGetStepDirection` functions to
change and query the direction of integration.

Added the function `MRIStepGetNumInnerStepperFails` to retrieve the number of
recoverable failures reported by the MRIStepInnerStepper.

Added a utility routine to wrap any valid ARKODE integrator for use as an
MRIStep inner stepper object, `ARKodeCreateMRIStepInnerStepper`.

The following DIRK schemes now have coefficients accurate to quad precision:

* `ARKODE_BILLINGTON_3_3_2`
* `ARKODE_KVAERNO_4_2_3`
* `ARKODE_CASH_5_2_4`
* `ARKODE_CASH_5_3_4`
* `ARKODE_KVAERNO_5_3_4`
* `ARKODE_KVAERNO_7_4_5`

#### CMake

The default value of `CMAKE_CUDA_ARCHITECTURES` is no longer set to `70` and is
now determined automatically by CMake. The previous default was only valid for
Volta GPUs while the automatically selected value will vary across compilers and
compiler versions. As such, users are encouraged to override this value with the
architecture for their system.

The build system has been updated to utilize the CMake LAPACK imported target
which should ease building SUNDIALS with LAPACK libraries that require setting
specific linker flags e.g., MKL.

#### Third Party Libraries

The Trilinos Teptra NVector interface has been updated to utilize CMake imported
targets added in Trilinos 14 to improve support for different Kokkos backends
with Trilinos. As such, Trilinos 14 or newer is required and the
`Trilinos_INTERFACE_*` CMake options have been removed.

Example programs using *hypre* have been updated to support v2.20 and newer.

### Bug Fixes

#### CMake

Fixed a CMake bug regarding usage of missing "print_warning" macro that was only
triggered when the deprecated `CUDA_ARCH` option was used.

Fixed a CMake configuration issue related to aliasing an `ALIAS` target when
using `ENABLE_KLU=ON` in combination with a static-only build of SuiteSparse.

Fixed a CMake issue which caused third-party CMake variables to be unset.  Users
may see more options in the CMake GUI now as a result of the fix.  See details
in GitHub Issue [#538](https://github.com/LLNL/sundials/issues/538).

#### NVector

Fixed a build failure with the SYCL NVector when using Intel oneAPI 2025.0
compilers. See GitHub Issue [#596](https://github.com/LLNL/sundials/issues/596).

Fixed compilation errors when building the Trilinos Teptra NVector with CUDA
support.

#### SUNMatrix

Fixed a [bug](https://github.com/LLNL/sundials/issues/581) in the sparse matrix
implementation of `SUNMatScaleAddI` which caused out of bounds writes unless
`indexvals` were in ascending order for each row/column.

#### SUNLinearSolver

Fixed a bug in the SPTFQMR linear solver where recoverable preconditioner errors
were reported as unrecoverable.

#### ARKODE

Fixed `ARKodeResize` not using the default `hscale` when an argument of `0` was
provided.

Fixed a memory leak that could occur if ``ARKodeSetDefaults`` is called
repeatedly.

Fixed the loading of ARKStep's default first order explicit method.

Fixed loading the default IMEX-MRI method if `ARKodeSetOrder` is used to specify
a third or fourth order method. Previously, the default second order method was
loaded in both cases.

Fixed potential memory leaks and out of bounds array accesses that could occur
in the ARKODE Lagrange interpolation module when changing the method order or
polynomial degree after re-initializing an integrator.

Fixed a bug in ARKODE when enabling rootfinding with fixed step sizes and the
initial value of the rootfinding function is zero. In this case, uninitialized
right-hand side data was used to compute a state value near the initial
condition to determine if any rootfinding functions are initially active.

Fixed a bug in MRIStep where the data supplied to the Hermite interpolation
module did not include contributions from the fast right-hand side
function. With this fix, users will see one additional fast right-hand side
function evaluation per slow step with the Hermite interpolation option.

Fixed a bug in SPRKStep when using compensated summations where the error vector
was not initialized to zero.

#### CVODE(S)

Fixed a bug where `CVodeSetProjFailEta` would ignore the `eta` parameter.

#### Fortran Interfaces

Fixed a bug in the 32-bit ``sunindextype`` Fortran interfaces to
``N_VGetSubvectorArrayPointer_ManyVector``,
``N_VGetSubvectorArrayPointer_MPIManyVector``, ``SUNBandMatrix_Column`` and
``SUNDenseMatrix_Column`` where 64-bit ``sunindextype`` interface functions were
used.

### Deprecation Notices

Deprecated the ARKStep-specific utility routine for wrapping an ARKStep instance
as an MRIStep inner stepper object, `ARKStepCreateMRIStepInnerStepper`. Use
`ARKodeCreateMRIStepInnerStepper` instead.

The ARKODE stepper specific functions to retrieve the number of right-hand side
function evaluations have been deprecated. Use `ARKodeGetNumRhsEvals` instead.

## Changes to SUNDIALS in release 7.1.1

### Bug Fixes

Fixed a [bug](https://github.com/LLNL/sundials/pull/523) in v7.1.0 with the SYCL
N_Vector `N_VSpace` function.

## Changes to SUNDIALS in release 7.1.0

### Major Features

Created shared user interface functions for ARKODE to allow more uniform control
over time-stepping algorithms, improved extensibility, and simplified code
maintenance. The corresponding stepper-specific user-callable functions are now
deprecated and will be removed in a future major release.

Added CMake infrastructure that enables externally maintained addons/plugins to
be *optionally* built with SUNDIALS. See the [Contributing
Guide](./CONTRIBUTING.md) for more details.

### New Features and Enhancements

Added support for Kokkos Kernels v4.

Added the following Runge-Kutta Butcher tables

* `ARKODE_FORWARD_EULER_1_1`
* `ARKODE_RALSTON_EULER_2_1_2`
* `ARKODE_EXPLICIT_MIDPOINT_EULER_2_1_2`
* `ARKODE_BACKWARD_EULER_1_1`
* `ARKODE_IMPLICIT_MIDPOINT_1_2`
* `ARKODE_IMPLICIT_TRAPEZOIDAL_2_2`

Added the following MRI coupling tables

* `ARKODE_MRI_GARK_FORWARD_EULER`
* `ARKODE_MRI_GARK_RALSTON2`
* `ARKODE_MRI_GARK_RALSTON3`
* `ARKODE_MRI_GARK_BACKWARD_EULER`
* `ARKODE_MRI_GARK_IMPLICIT_MIDPOINT`
* `ARKODE_IMEX_MRI_GARK_EULER`
* `ARKODE_IMEX_MRI_GARK_TRAPEZOIDAL`
* `ARKODE_IMEX_MRI_GARK_MIDPOINT`

Added `ARKodeButcherTable_ERKIDToName` and `ARKodeButcherTable_DIRKIDToName` to
convert a Butcher table ID to a string representation.

Added the function ``ARKodeSetAutonomous`` in ARKODE to indicate that the
implicit right-hand side function does not explicitly depend on time. When using
the trivial predictor, an autonomous problem may reuse implicit function
evaluations across stage solves to reduce the total number of function
evaluations.

Users may now disable interpolated output in ARKODE by passing `ARK_INTERP_NONE`
to `ARKodeSetInterpolantType`. When interpolation is disabled, rootfinding is
not supported, implicit methods must use the trivial predictor (the default
option), and interpolation at stop times cannot be used (interpolating at stop
times is disabled by default). With interpolation disabled, calling
`ARKodeEvolve` in `ARK_NORMAL` mode will return at or past the requested output
time (setting a stop time may still be used to halt the integrator at a specific
time). Disabling interpolation will reduce the memory footprint of an integrator
by two or more state vectors (depending on the interpolant type and degree)
which can be beneficial when interpolation is not needed e.g., when integrating
to a final time without output in between or using an explicit fast time scale
integrator with an MRI method.

Added "Resize" capability to ARKODE's SPRKStep time-stepping module.

Enabled the Fortran interfaces to build with 32-bit `sunindextype`.

### Bug Fixes

Updated the CMake variable `HIP_PLATFORM` default to `amd` as the previous
default, `hcc`, is no longer recognized in ROCm 5.7.0 or newer. The new default
is also valid in older version of ROCm (at least back to version 4.3.1).

Renamed the DPCPP value for the `SUNDIALS_GINKGO_BACKENDS` CMake option to `SYCL`
to match Ginkgo's updated naming convention.

Changed the CMake version compatibility mode for SUNDIALS to `AnyNewerVersion`
instead of `SameMajorVersion`. This fixes the issue seen
[here](https://github.com/AMReX-Codes/amrex/pull/3835).

Fixed a CMake bug that caused an MPI linking error for our C++ examples in some
instances. Fixes [GitHub Issue #464](https://github.com/LLNL/sundials/issues/464).

Fixed the runtime library installation path for windows systems. This fix
changes the default library installation path from
`CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_LIBDIR` to
`CMAKE_INSTALL_PREFIX/CMAKE_INSTALL_BINDIR`.

Fixed conflicting `.lib` files between shared and static libs when using `MSVC`
on Windows.

Fixed invalid `SUNDIALS_EXPORT` generated macro when building both shared and
static libs.

Fixed a bug in some Fortran examples where `c_null_ptr` was passed as an
argument to a function pointer instead of `c_null_funptr`. This caused
compilation issues with the Cray Fortran compiler.

Fixed a bug in the HIP execution policies where `WARP_SIZE` would not be set
with ROCm 6.0.0 or newer.

Fixed a bug that caused error messages to be cut off in some cases. Fixes
[GitHub Issue #461](https://github.com/LLNL/sundials/issues/461).

Fixed a memory leak when an error handler was added to a `SUNContext`. Fixes
[GitHub Issue #466](https://github.com/LLNL/sundials/issues/466).

Fixed a bug where `MRIStepEvolve` would not handle a recoverable error produced
from evolving the inner stepper.

Added missing `SetRootDirection` and `SetNoInactiveRootWarn` functions to
ARKODE's SPRKStep time-stepping module.

Fixed a bug in `ARKodeSPRKTable_Create` where the coefficient arrays were not
allocated.

Fix bug on LLP64 platforms (like Windows 64-bit) where `KLU_INDEXTYPE` could be
32 bits wide even if `SUNDIALS_INT64_T` is defined.

Check if size of `SuiteSparse_long` is 8 if the size of `sunindextype` is 8
when using KLU.

Fixed several build errors with the Fortran interfaces on Windows systems.

### Deprecation Notices

Numerous ARKODE stepper-specific functions are now deprecated in favor of
ARKODE-wide functions.

Deprecated the `ARKStepSetOptimalParams` function. Since this function does not have an
ARKODE-wide equivalent, instructions have been added to the user guide for how
to retain the current functionality using other user-callable functions.

The unsupported implementations of `N_VGetArrayPointer` and `N_VSetArrayPointer`
for the *hypre* and PETSc vectors are now deprecated. Users should access the
underlying wrapped external library vector objects instead with
`N_VGetVector_ParHyp` and `N_VGetVector_Petsc`, respectively.

## Changes to SUNDIALS in release v7.0.0

### Major Feature

SUNDIALS now has more robust and uniform error handling. Non-release builds will
be built with additional error checking by default. See the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
section in the user guide for details.

### Breaking Changes

#### Minimum C Standard

SUNDIALS now requires using a compiler that supports a subset of the C99
standard. Note with the Microsoft C/C++ compiler the subset of C99 features
utilized by SUNDIALS are available starting with [Visual Studio 2015](https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170#c-standard-library-features-1).

#### Minimum CMake Version

CMake 3.18 or newer is now required when building SUNDIALS.

#### Deprecated Types and Functions Removed

The previously deprecated types `realtype` and `booleantype` were removed from
`sundials_types.h` and replaced with `sunrealtype` and `sunbooleantype`. The
deprecated names for these types can be used by including the header file
`sundials_types_deprecated.h` but will be fully removed in the next major
release. Functions, types and header files that were previously deprecated have
also been removed.

#### Error Handling Changes

With the addition of the new error handling capability, the `*SetErrHandlerFn`
and `*SetErrFile` functions in CVODE(S), IDA(S), ARKODE, and KINSOL have been
removed. Users of these functions can use the functions
`SUNContext_PushErrHandler`, and `SUNLogger_SetErrorFilename` instead. For
further details see the
[Error Checking](https://sundials.readthedocs.io/en/latest/sundials/Errors_link.html)
and
[Logging](https://sundials.readthedocs.io/en/latest/sundials/Logging_link.html)
sections in the documentation.

In addition the following names/symbols were replaced by `SUN_ERR_*` codes:

| Removed                        | Replaced with `SUNErrCode`        |
|:-------------------------------|:----------------------------------|
| `SUNLS_SUCCESS`                | `SUN_SUCCESS`                     |
| `SUNLS_UNRECOV_FAILURE`        | no replacement (value was unused) |
| `SUNLS_MEM_NULL`               | `SUN_ERR_ARG_CORRUPT`             |
| `SUNLS_ILL_INPUT`              | `SUN_ERR_ARG_*`                   |
| `SUNLS_MEM_FAIL`               | `SUN_ERR_MEM_FAIL`                |
| `SUNLS_PACKAGE_FAIL_UNREC`     | `SUN_ERR_EXT_FAIL`                |
| `SUNLS_VECTOROP_ERR`           | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_SUCCESS`              | `SUN_SUCCESS`                     |
| `SUN_NLS_MEM_NULL`             | `SUN_ERR_ARG_CORRUPT`             |
| `SUN_NLS_MEM_FAIL`             | `SUN_ERR_MEM_FAIL`                |
| `SUN_NLS_ILL_INPUT`            | `SUN_ERR_ARG_*`                   |
| `SUN_NLS_VECTOROP_ERR`         | `SUN_ERR_OP_FAIL`                 |
| `SUN_NLS_EXT_FAIL`             | `SUN_ERR_EXT_FAIL`                |
| `SUNMAT_SUCCESS`               | `SUN_SUCCESS`                     |
| `SUNMAT_ILL_INPUT`             | `SUN_ERR_ARG_*`                   |
| `SUNMAT_MEM_FAIL`              | `SUN_ERR_MEM_FAIL`                |
| `SUNMAT_OPERATION_FAIL`        | `SUN_ERR_OP_FAIL`                 |
| `SUNMAT_MATVEC_SETUP_REQUIRED` | `SUN_ERR_OP_FAIL`                 |

The following functions have had their signature updated to ensure they can
leverage the new SUNDIALS error handling capabilities.

```c
// From sundials_futils.h
SUNDIALSFileOpen
SUNDIALSFileClose

// From sundials_memory.h
SUNMemoryNewEmpty
SUNMemoryHelper_Alias
SUNMemoryHelper_Wrap

// From sundials_nvector.h
N_VNewVectorArray
```

#### SUNComm Type Added

We have replaced the use of a type-erased (i.e., `void*`) pointer to a
communicator in place of `MPI_Comm` throughout the SUNDIALS API with a
`SUNComm`, which is just a typedef to an `int` in builds without MPI
and a typedef to a `MPI_Comm` in builds with MPI. As a result:

* When MPI is enabled, all SUNDIALS libraries will include MPI symbols and
  applications will need to include the path for MPI headers and link against
  the corresponding MPI library.

* All users will need to update their codes because the call to
  `SUNContext_Create` now takes a `SUNComm` instead
  of type-erased pointer to a communicator. For non-MPI codes,
  pass `SUN_COMM_NULL` to the `comm` argument instead of
  `NULL`. For MPI codes, pass the `MPI_Comm` directly.

* The same change must be made for calls to
  `SUNLogger_Create` or `SUNProfiler_Create`.

* Some users will need to update their calls to `N_VGetCommunicator`, and
  update any custom `N_Vector` implementations that provide
  `N_VGetCommunicator`, since it now returns a `SUNComm`.

The change away from type-erased pointers for `SUNComm` fixes problems like the
one described in [GitHub Issue #275](https://github.com/LLNL/sundials/issues/275).

The SUNLogger is now always MPI-aware if MPI is enabled in SUNDIALS and the
`SUNDIALS_LOGGING_ENABLE_MPI` CMake option and macro definition were removed
accordingly.

#### SUNDIALS Core Library

Users now need to link to `sundials_core` in addition to the libraries already
linked to. This will be picked up automatically in projects that use the
SUNDIALS CMake target. The library `sundials_generic` has been superseded by
`sundials_core` and is no longer available. This fixes some duplicate symbol
errors on Windows when linking to multiple SUNDIALS libraries.

#### Fortran Interface Modules Streamlined

We have streamlined the Fortran modules that need to be included by users by combining
the SUNDIALS core into one Fortran module, `fsundials_core_mod`. Modules for
implementations of the core APIs still exist (e.g., for the Dense linear solver there
is `fsunlinsol_dense_mod`) as do the modules for the SUNDIALS packages (e.g., `fcvode_mod`).
The following modules are the ones that have been consolidated into `fsundials_core_mod`:

```
fsundials_adaptcontroller_mod
fsundials_context_mod
fsundials_futils_mod
fsundials_linearsolver_mod
fsundials_logger_mod
fsundials_matrix_mod
fsundials_nonlinearsolver_mod
fsundials_nvector_mod
fsundials_profiler_mod
fsundials_types_mod
```

### Minor Changes

The `CMAKE_BUILD_TYPE` defaults to `RelWithDebInfo` mode now i.e., SUNDIALS
will be built with optimizations and debugging symbols enabled by default.
Previously the build type was unset by default so no optimization or debugging
flags were set.

The advanced CMake options to override the inferred LAPACK name-mangling scheme
have been updated from `SUNDIALS_F77_FUNC_CASE` and
`SUNDIALS_F77_FUNC_UNDERSCORES` to `SUNDIALS_LAPACK_CASE` and
`SUNDIALS_LAPACK_UNDERSCORES`, respectively.

As a subset of C99 is now required the CMake option `USE_GENERIC_MATH` as been
removed.

The C++ convenience classes (e.g., `sundials::Context`) have been moved to
from SUNDIALS `.h` headers to corresponding `.hpp` headers (e.g.,
`sundials/sundials_context.hpp`) so C++ codes do not need to compile with
C++14 support when using the C API.

Converted most previous Fortran 77 and 90 examples to use SUNDIALS' Fortran 2003
interface.

### Bug Fixes

Fixed [#329](https://github.com/LLNL/sundials/issues/329) so that C++20 aggregate initialization can be used.

Fixed integer overflow in the internal SUNDIALS hashmap. This resolves
[#409](https://github.com/LLNL/sundials/issues/409) and
[#249](https://github.com/LLNL/sundials/issues/249).

### Deprecation Notice

The functions in `sundials_math.h` will be deprecated in the next release.

```c
  sunrealtype SUNRpowerI(sunrealtype base, int exponent);
  sunrealtype SUNRpowerR(sunrealtype base, sunrealtype exponent);
  sunbooleantype SUNRCompare(sunrealtype a, sunrealtype b);
  sunbooleantype SUNRCompareTol(sunrealtype a, sunrealtype b, sunrealtype tol);
  sunrealtype SUNStrToReal(const char* str);
```

Additionally, the following header files (and everything in them) will be
deprecated -- users who rely on these are recommended to transition to the
corresponding `SUNMatrix` and `SUNLinearSolver` modules:

```c
sundials_direct.h
sundials_dense.h
sundials_band.h
```

## Changes to SUNDIALS in release 6.7.0

### Major Feature

Added the `SUNAdaptController` base class, ported ARKODE's internal
implementations of time step controllers into implementations of this class,
and updated ARKODE to use these objects instead of its own implementations.
Added `ARKStepSetAdaptController` and `ERKStepSetAdaptController` routines
so that users can modify controller parameters, or even provide custom
implementations.

### New Features

Improved computational complexity of `SUNMatScaleAddI_Sparse` from `O(M*N)` to
`O(NNZ)`.

Added Fortran support for the LAPACK dense `SUNLinearSolver` implementation.

Added the routines `ARKStepSetAdaptivityAdjustment` and
`ERKStepSetAdaptivityAdjustment`, that allow users to adjust the
value for the method order supplied to the temporal adaptivity controllers.
The ARKODE default for this adjustment has been -1 since its initial
release, but for some applications a value of 0 is more appropriate.
Users who notice that their simulations encounter a large number of
temporal error test failures may want to experiment with adjusting this value.

Added the third order ERK method `ARKODE_SHU_OSHER_3_2_3`, the fourth order
ERK method `ARKODE_SOFRONIOU_SPALETTA_5_3_4`, the sixth order ERK method
`ARKODE_VERNER_9_5_6`, the seventh order ERK method `ARKODE_VERNER_10_6_7`,
the eighth order ERK method `ARKODE_VERNER_13_7_8`, and the ninth order ERK
method `ARKODE_VERNER_16_8_9`.

ARKStep, ERKStep, MRIStep, and SPRKStep were updated to remove a potentially
unnecessary right-hand side evaluation at the end of an integration. ARKStep was
additionally updated to remove extra right-hand side evaluations when using an
explicit method or an implicit method with an explicit first stage.

The `MRIStepInnerStepper` class in MRIStep was updated to make supplying an
`MRIStepInnerFullRhsFn` optional.

### Bug Fixes

Changed the `SUNProfiler` so that it does not rely on `MPI_WTime` in any case.
This fixes [GitHub Issue #312](https://github.com/LLNL/sundials/issues/312).

Fixed scaling bug in `SUNMatScaleAddI_Sparse` for non-square matrices.

Fixed a regression introduced by the stop time bug fix in v6.6.1 where ARKODE,
CVODE, CVODES, IDA, and IDAS would return at the stop time rather than the
requested output time if the stop time was reached in the same step in which the
output time was passed.

Fixed a bug in ERKStep where methods with `c[s-1] = 1` but `a[s-1,j] != b[j]`
were incorrectly treated as having the first same as last (FSAL) property.

Fixed a bug in ARKODE where `ARKStepSetInterpolateStopTime` would return an
interpolated solution at the stop time in some cases when interpolation was
disabled.

Fixed a bug in `ARKStepSetTableNum` wherein it did not recognize
`ARKODE_ARK2_ERK_3_1_2` and `ARKODE_ARK2_DIRK_3_1_2` as a valid additive
Runge-Kutta Butcher table pair.

Fixed a bug in `MRIStepCoupling_Write` where explicit coupling tables were not
written to the output file pointer.

Fixed missing soversions in some `SUNLinearSolver` and `SUNNonlinearSolver`
CMake targets.

Renamed some internal types in CVODES and IDAS to allow both packages to be
built together in the same binary.

## Changes to SUNDIALS in release 6.6.2

Fixed the build system support for MAGMA when using a NVIDIA HPC SDK
installation of CUDA and fixed the targets used for rocBLAS and rocSPARSE.

## Changes to SUNDIALS in release 6.6.1

### New Features

Updated the Tpetra NVector interface to support Trilinos 14.

### Bug Fixes

Fixed a memory leak when destroying a CUDA, HIP, SYCL, or system SUNMemoryHelper
object.

Fixed a bug in ARKODE, CVODE, CVODES, IDA, and IDAS where the stop time may not
be cleared when using normal mode if the requested output time is the same as
the stop time. Additionally, with ARKODE, CVODE, and CVODES this fix removes an
unnecessary interpolation of the solution at the stop time that could occur in
this case.

## Changes to SUNDIALS in release 6.6.0

### Major Features

A new time-stepping module, `SPRKStep`, was added to ARKODE. This time-stepper
provides explicit symplectic partitioned Runge-Kutta methods up to order 10
for separable Hamiltonian systems.

Added support for relaxation Runge-Kutta methods in ERKStep and ARKStep in
ARKODE.

### New Features

Updated CVODE, CVODES and ARKODE default behavior when returning the solution when
the internal time has reached a user-specified stop time. Previously, the output
solution was interpolated to the value of `tstop`; the default is now to copy the
internal solution vector. Users who wish to revert to interpolation may call a new
routine `CVodeSetInterpolateStopTime`, `ARKStepSetInterpolateStopTime`,
`ERKStepSetInterpolateStopTime`, or `MRIStepSetInterpolateStopTime`.

Added the second order IMEX method from Giraldo, Kelly, and Constantinescu 2013
as the default second order IMEX method in ARKStep. The explicit table is given
by `ARKODE_ARK2_ERK_3_1_2` and the implicit table by `ARKODE_ARK2_DIRK_3_1_2`.

Updated the F2003 utility routines `SUNDIALSFileOpen` and `SUNDIALSFileClose`
to support user specification of `stdout` and `stderr` strings for the output
file names.

## Bug Fixes

A potential bug was fixed when using inequality constraint handling and
calling `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` after a failed
step in which an inequality constraint violation occurred. In this case, the
values returned by `ARKStepGetEstLocalErrors` or `ERKStepGetEstLocalErrors` may
have been invalid.

## Changes to SUNDIALS in release 6.5.1

### New Features

Added the functions `ARKStepClearStopTime`, `ERKStepClearStopTime`,
`MRIStepClearStopTime`, `CVodeClearStopTime`, and `IDAClearStopTime` to
disable a previously set stop time.

The default interpolant in ARKODE when using a first order method has been
updated to a linear interpolant to ensure values obtained by the integrator are
returned at the ends of the time interval. To restore the previous behavior of
using a constant interpolant call `ARKStepSetInterpolantDegree`,
`ERKStepSetInterpolantDegree`, or `MRIStepSetInterpolantDegree` and set the
interpolant degree to zero before evolving the problem.

### Bug Fixes

Fixed build errors when using SuperLU_DIST with ROCM enabled to target AMD GPUs.

Fixed compilation errors in some SYCL examples when using the `icx` compiler.

## Changes to SUNDIALS in release 6.5.0

### New Features

A new capability to keep track of memory allocations made through the `SUNMemoryHelper`
classes has been added. Memory allocation stats can be accessed through the
`SUNMemoryHelper_GetAllocStats` function. See the documentation for
the `SUNMemoryHelper` classes for more details.

Added the functions `ARKStepGetJac`, `ARKStepGetJacTime`,
`ARKStepGetJacNumSteps`, `MRIStepGetJac`, `MRIStepGetJacTime`,
`MRIStepGetJacNumSteps`, `CVodeGetJac`, `CVodeGetJacTime`,
`CVodeGetJacNumSteps`, `IDAGetJac`, `IDAGetJacCj`, `IDAGetJacTime`,
`IDAGetJacNumSteps`, `KINGetJac`, `KINGetJacNumIters` to assist in
debugging simulations utilizing a matrix-based linear solver.

Added support for CUDA 12.

Added support for the SYCL backend with RAJA 2022.x.y.

### Bug Fixes

Fixed an underflow bug during root finding in ARKODE, CVODE, CVODES, IDA and
IDAS. This fixes [GitHub Issue #57](https://github.com/LLNL/sundials/issues/57>).

Fixed an issue with finding oneMKL when using the `icpx` compiler with the
`-fsycl` flag as the C++ compiler instead of `dpcpp`.

Fixed the shape of the arrays returned by `FN_VGetArrayPointer` functions as well
as the `FSUNDenseMatrix_Data`, `FSUNBandMatrix_Data`, `FSUNSparseMatrix_Data`,
`FSUNSparseMatrix_IndexValues`, and `FSUNSparseMatrix_IndexPointers` functions.
Compiling and running code that uses the SUNDIALS Fortran interfaces with
bounds checking will now work.

Fixed an implicit conversion error in the Butcher table for ESDIRK5(4)7L[2]SA2.

## Changes to SUNDIALS in release 6.4.1

Fixed a bug with the Kokkos interfaces that would arise when using clang.

Fixed a compilation error with the Intel oneAPI 2022.2 Fortran compiler in the
Fortran 2003 interface test for the serial `N_Vector`.

Fixed a bug in the SUNLINSOL_LAPACKBAND and SUNLINSOL_LAPACKDENSE modules
which would cause the tests to fail on some platforms.

## Changes to SUNDIALS in release 6.4.0

### New Requirements

CMake 3.18.0 or newer is now required for CUDA support.

A C++14 compliant compiler is now required for C++ based features and examples
e.g., CUDA, HIP, RAJA, Trilinos, SuperLU_DIST, MAGMA, Ginkgo, and Kokkos.

### Major Features

Added support for the [Ginkgo](https://ginkgo-project.github.io/) linear algebra
library. This support includes new `SUNMatrix` and `SUNLinearSolver`
implementations, see the `SUNMATRIX_GINKGO` and `SUNLINEARSOLVER_GINKGO`
sections in the documentation for more information.

Added new `NVector`, dense `SUNMatrix`, and dense `SUNLinearSolver`
implementations utilizing [Kokkos Ecosystem](https://kokkos.org/) for
performance portability, see the `NVECTOR_KOKKOS`, `SUNMATRIX_KOKKOSDENSE` and
`SUNLINEARSOLVER_KOKKOSDENSE` sections in the documentation for more
information.

### New Features

Added support for GPU enabled SuperLU_DIST and SuperLU_DIST v8.x.x. Removed
support for SuperLU_DIST v6.x.x or older. Fix mismatched definition and
declaration bug in SuperLU_DIST matrix constructor.

Added the functions `ARKStepSetTableName`, `ERKStepSetTableName`,
`MRIStepCoupling_LoadTableByName`, `ARKodeButcherTable_LoadDIRKByName`, and
`ARKodeButcherTable_LoadERKByName` to load a table from a string.

### Bug Fixes

Fixed a bug in the CUDA and HIP vectors where `N_VMaxNorm` would return the
minimum positive floating-point value for the zero vector.

Fixed memory leaks/out of bounds memory accesses in the ARKODE MRIStep module
that could occur when attaching a coupling table after reinitialization with a
different number of stages than originally selected.

Fixed a memory leak in CVODE and CVODES where the projection memory would not be
deallocated when calling `CVodeFree`.

## Changes to SUNDIALS in release 6.3.0

### New Features

Added `GetUserData` functions in each package to retrieve the user data pointer
provided to `SetUserData` functions. See `ARKStepGetUserData`,
`ERKStepGetUserData`, `MRIStepGetUserData`, `CVodeGetUserData`,
`IDAGetUserData`, or `KINGetUserData` for more information.

Added a variety of embedded DIRK methods from [Kennedy & Carpenter,
NASA TM-2016-219173, 2016] and [Kennedy & Carpenter, Appl. Numer. Math., 146, 2019] to
ARKODE.

Updated `MRIStepReset` to call the corresponding `MRIStepInnerResetFn` with the same
(*tR*,*yR*) arguments for the `MRIStepInnerStepper` object that is used to evolve the
MRI "fast" time scale subproblems.

Added a new [example](examples/cvode/serial/cvRocket_dns.c) which
demonstrates using CVODE with a discontinuous right-hand-side function
and rootfinding.

### Bug Fixes

Fixed a bug in `ERKStepReset`, `ERKStepReInit`, `ARKStepReset`, `ARKStepReInit`,
`MRIStepReset`, and `MRIStepReInit` where a previously-set value of *tstop* (from
a call to `ERKStepSetStopTime`, `ARKStepSetStopTime`, or `MRIStepSetStopTime`,
respectively) would not be cleared.

Fixed the unituitive behavior of the `USE_GENERIC_MATH` CMake option which
caused the double precision math functions to be used regardless of the value of
`SUNDIALS_PRECISION`. Now, SUNDIALS will use precision appropriate math
functions when they are available and the user may provide the math library to
link to via the advanced CMake option `SUNDIALS_MATH_LIBRARY`.

Changed `SUNDIALS_LOGGING_ENABLE_MPI` CMake option default to be 'OFF'. This
fixes [GitHub Issue #177](https://github.com/LLNL/sundials/issues/177).

## Changes to SUNDIALS in release 6.2.0

### Major Features

Added the `SUNLogger` API which provides a SUNDIALS-wide
mechanism for logging of errors, warnings, informational output,
and debugging output.

Added support to CVODES for integrating IVPs with constraints using BDF methods
and projecting the solution onto the constraint manifold with a user defined
projection function. This implementation is accompanied by additions to the
CVODES user documentation and examples.

### New Features

Added the function `SUNProfiler_Reset` to reset the region timings and counters
to zero.

Added the following functions to output all of the integrator, nonlinear solver,
linear solver, and other statistics in one call:

* `ARKStepPrintAllStats`
* `ERKStepPrintAllStats`
* `MRIStepPrintAllStats`
* `CVodePrintAllStats`
* `IDAPrintAllStats`
* `KINPrintAllStats`

The file `scripts/sundials_csv.py` contains functions for parsing the
comma-separated value (CSV) output files when using the CSV output format.

Added functions to CVODE, CVODES, IDA, and IDAS to change the default step size
adaptivity parameters. For more information see the documentation for:

* `CVodeSetEtaFixedStepBounds`
* `CVodeSetEtaMaxFirstStep`
* `CVodeSetEtaMaxEarlyStep`
* `CVodeSetNumStepsEtaMaxEarlyStep`
* `CVodeSetEtaMax`
* `CVodeSetEtaMin`
* `CVodeSetEtaMinErrFail`
* `CVodeSetEtaMaxErrFail`
* `CVodeSetNumFailsEtaMaxErrFail`
* `CVodeSetEtaConvFail`
* `IDASetEtaFixedStepBounds`
* `IDAsetEtaMax`
* `IDASetEtaMin`
* `IDASetEtaLow`
* `IDASetEtaMinErrFail`
* `IDASetEtaConvFail`

Added the functions `ARKStepSetDeduceImplicitRhs` and
`MRIStepSetDeduceImplicitRhs` to optionally remove an evaluation of the implicit
right-hand side function after nonlinear solves. See the mathematical
considerations section of the user guide for information on using this
optimization.

Added the function `MRIStepSetOrder` to select the default MRI method of a given
order.

Added the functions `CVodeSetDeltaGammaMaxLSetup` and
`CVodeSetDeltaGammaMaxBadJac` in CVODE and CVODES to adjust the `gamma` change
thresholds to require a linear solver setup or Jacobian/precondition update,
respectively.

Added the function `IDASetDetlaCjLSetup` in IDA and IDAS to adjust the parameter
that determines when a change in `c_j` requires calling the linear solver setup
function.

Added the function `IDASetMinStep` to set a minimum step size.

### Bug Fixes

Fixed the `SUNContext` convenience class for C++ users to disallow copy
construction and allow move construction.

The behavior of `N_VSetKernelExecPolicy_Sycl` has been updated to be consistent
with the CUDA and HIP vectors. The input execution policies are now cloned and
may be freed after calling `N_VSetKernelExecPolicy_Sycl`. Additionally, `NULL`
inputs are now allowed and, if provided, will reset the vector execution
policies to the defaults.

A memory leak in the SYCL vector was fixed where the execution policies were not
freed when the vector was destroyed.

The include guard in `nvector_mpimanyvector.h` has been corrected to enable
using both the ManyVector and MPIManyVector vector implementations in the same
simulation.

A bug was fixed in the ARKODE, CVODE(S), and IDA(S) functions to retrieve the
number of nonlinear solver failures. The failure count returned was the number
of failed *steps* due to a nonlinear solver failure i.e., if a nonlinear solve
failed with a stale Jacobian or preconditioner but succeeded after updating the
Jacobian or preconditioner, the initial failure was not included in the
nonlinear solver failure count. The following functions have been updated to
return the total number of nonlinear solver failures:

* `ARKStepGetNumNonlinSolvConvFails`
* `ARKStepGetNonlinSolvStats`
* `MRIStepGetNumNonlinSolvConvFails`
* `MRIStepGetNonlinSolvStats`
* `CVodeGetNumNonlinSolvConvFails`
* `CVodeGetNonlinSolvStats`
* `CVodeGetSensNumNonlinSolvConvFails`
* `CVodeGetSensNonlinSolvStats`
* `CVodeGetStgrSensNumNonlinSolvConvFails`
* `CVodeGetStgrSensNonlinSolvStats`
* `IDAGetNumNonlinSolvConvFails`
* `IDAGetNonlinSolvStats`
* `IDAGetSensNumNonlinSolvConvFails`
* `IDAGetSensNonlinSolvStats`

As a result of this change users may see an increase in the number of failures
reported from the above functions. The following functions have been added to
retrieve the number of failed steps due to a nonlinear solver failure i.e., the
counts previously returned by the above functions:

* `ARKStepGetNumStepSolveFails`
* `MRIStepGetNumStepSolveFails`
* `CVodeGetNumStepSolveFails`
* `CVodeGetNumStepSensSolveFails`
* `CVodeGetNumStepStgrSensSolveFails`
* `IDAGetNumStepSolveFails`
* `IDAGetNumStepSensSolveFails`

Changed exported SUNDIALS PETSc CMake targets to be INTERFACE IMPORTED instead
of UNKNOWN IMPORTED.

### Deprecation Notice

Deprecated the following functions, it is recommended to use the `SUNLogger` API
instead.

* `ARKStepSetDiagnostics`
* `ERKStepSetDiagnostics`
* `MRIStepSetDiagnostics`
* `KINSetInfoFile`
* `SUNNonlinSolSetPrintLevel_Newton`
* `SUNNonlinSolSetInfoFile_Newton`
* `SUNNonlinSolSetPrintLevel_FixedPoint`
* `SUNNonlinSolSetInfoFile_FixedPoint`
* `SUNLinSolSetInfoFile_PCG`
* `SUNLinSolSetPrintLevel_PCG`
* `SUNLinSolSetInfoFile_SPGMR`
* `SUNLinSolSetPrintLevel_SPGMR`
* `SUNLinSolSetInfoFile_SPFGMR`
* `SUNLinSolSetPrintLevel_SPFGMR`
* `SUNLinSolSetInfoFile_SPTFQM`
* `SUNLinSolSetPrintLevel_SPTFQMR`
* `SUNLinSolSetInfoFile_SPBCGS`
* `SUNLinSolSetPrintLevel_SPBCGS`

The `SUNLinSolSetInfoFile_*` and `SUNNonlinSolSetInfoFile_*` family of
functions are now enabled by setting the CMake option `SUNDIALS_LOGGING_LEVEL`
to a value `>= 3`.

## Changes to SUNDIALS in release 6.1.1

### New Features

Added new Fortran example program,
`examples/arkode/F2003_serial/ark_kpr_mri_f2003.f90` demonstrating MRI
capabilities.

### Bug Fixes

Fixed exported `SUNDIALSConfig.cmake`.

Fixed Fortran interface to `MRIStepInnerStepper` and `MRIStepCoupling`
structures and functions.

## Changes to SUNDIALS in release 6.1.0

### New Features

Added new reduction implementations for the CUDA and HIP vectors that use
shared memory (local data storage) instead of atomics. These new implementations
are recommended when the target hardware does not provide atomic support for the
floating point precision that SUNDIALS is being built with. The HIP vector uses
these by default, but the `N_VSetKernelExecPolicy_Cuda` and
`N_VSetKernelExecPolicy_Hip` functions can be used to choose between
different reduction implementations.

`SUNDIALS::<lib>` targets with no static/shared suffix have been added for use
within the build directory (this mirrors the targets exported on installation).

`CMAKE_C_STANDARD` is now set to 99 by default.

### Bug Fixes

Fixed exported `SUNDIALSConfig.cmake` when profiling is enabled without Caliper.

Fixed `sundials_export.h` include in `sundials_config.h`.

Fixed memory leaks in the SuperLU_MT linear solver interface.

## Changes to SUNDIALS in release 6.0.0

### Breaking Changes

#### SUNContext Object Added

SUNDIALS v6.0.0 introduces a new `SUNContext` object on which all other SUNDIALS
objects depend. As such, the constructors for all SUNDIALS packages, vectors,
matrices, linear solvers, nonlinear solvers, and memory helpers have been
updated to accept a context as the last input. Users upgrading to SUNDIALS
v6.0.0 will need to call `SUNContext_Create` to create a context object with
before calling any other SUNDIALS library function, and then provide this object
to other SUNDIALS constructors. The context object has been introduced to allow
SUNDIALS to provide new features, such as the profiling/instrumentation also
introduced in this release, while maintaining thread-safety. See the
documentation section on the `SUNContext` for more details.

The script `upgrade-to-sundials-6-from-5.sh` has been provided with this release
(and obtainable from the GitHub release page) to help ease the transition to
SUNDIALS v6.0.0. The script will add a `SUNCTX_PLACEHOLDER` argument to all of
the calls to SUNDIALS constructors that now require a `SUNContext` object. It
can also update deprecated SUNDIALS constants/types to the new names. It can be
run like this:

```
> ./upgrade-to-sundials-6-from-5.sh <files to update>
```

#### Updated SUNMemoryHelper Function Signatures

The `SUNMemoryHelper` functions `Alloc`, `Dealloc`, and `Copy` have been updated
to accept an opaque handle as the last input. At a minimum, existing
`SUNMemoryHelper` implementations will need to update these functions to accept
the additional argument. Typically, this handle is the execution stream (e.g., a
CUDA/HIP stream or SYCL queue) for the operation. The CUDA, HIP, and SYCL
`SUNMemoryHelper` implementations have been updated accordingly. Additionally,
the constructor for the SYCL implementation has been updated to remove the SYCL
queue as an input.

#### Deprecated Functions Removed

The previously deprecated constructor `N_VMakeWithManagedAllocator_Cuda` and
the function `N_VSetCudaStream_Cuda` have been removed and replaced with
`N_VNewWithMemHelp_Cuda` and `N_VSetKernelExecPolicy_Cuda` respectively.

The previously deprecated macros `PVEC_REAL_MPI_TYPE` and
`PVEC_INTEGER_MPI_TYPE` have been removed and replaced with
`MPI_SUNREALTYPE` and `MPI_SUNINDEXTYPE` respectively.

The following previously deprecated functions have been removed

| Removed                   | Replaced with                    |
|:--------------------------|:---------------------------------|
| `SUNBandLinearSolver`     | `SUNLinSol_Band`                 |
| `SUNDenseLinearSolver`    | `SUNLinSol_Dense`                |
| `SUNKLU`                  | `SUNLinSol_KLU`                  |
| `SUNKLUReInit`            | `SUNLinSol_KLUReInit`            |
| `SUNKLUSetOrdering`       | `SUNLinSol_KLUSetOrdering`       |
| `SUNLapackBand`           | `SUNLinSol_LapackBand`           |
| `SUNLapackDense`          | `SUNLinSol_LapackDense`          |
| `SUNPCG`                  | `SUNLinSol_PCG`                  |
| `SUNPCGSetPrecType`       | `SUNLinSol_PCGSetPrecType`       |
| `SUNPCGSetMaxl`           | `SUNLinSol_PCGSetMaxl`           |
| `SUNSPBCGS`               | `SUNLinSol_SPBCGS`               |
| `SUNSPBCGSSetPrecType`    | `SUNLinSol_SPBCGSSetPrecType`    |
| `SUNSPBCGSSetMaxl`        | `SUNLinSol_SPBCGSSetMaxl`        |
| `SUNSPFGMR`               | `SUNLinSol_SPFGMR`               |
| `SUNSPFGMRSetPrecType`    | `SUNLinSol_SPFGMRSetPrecType`    |
| `SUNSPFGMRSetGSType`      | `SUNLinSol_SPFGMRSetGSType`      |
| `SUNSPFGMRSetMaxRestarts` | `SUNLinSol_SPFGMRSetMaxRestarts` |
| `SUNSPGMR`                | `SUNLinSol_SPGMR`                |
| `SUNSPGMRSetPrecType`     | `SUNLinSol_SPGMRSetPrecType`     |
| `SUNSPGMRSetGSType`       | `SUNLinSol_SPGMRSetGSType`       |
| `SUNSPGMRSetMaxRestarts`  | `SUNLinSol_SPGMRSetMaxRestarts`  |
| `SUNSPTFQMR`              | `SUNLinSol_SPTFQMR`              |
| `SUNSPTFQMRSetPrecType`   | `SUNLinSol_SPTFQMRSetPrecType`   |
| `SUNSPTFQMRSetMaxl`       | `SUNLinSol_SPTFQMRSetMaxl`       |
| `SUNSuperLUMT`            | `SUNLinSol_SuperLUMT`            |
| `SUNSuperLUMTSetOrdering` | `SUNLinSol_SuperLUMTSetOrdering` |

The deprecated functions `MRIStepGetCurrentButcherTables` and
`MRIStepWriteButcher` and the utility functions `MRIStepSetTable` and
`MRIStepSetTableNum` have been removed. Users wishing to create an MRI-GARK
method from a Butcher table should use `MRIStepCoupling_MIStoMRI` to create
the corresponding MRI coupling table and attach it with `MRIStepSetCoupling`.

The previously deprecated functions `ARKStepSetMaxStepsBetweenLSet` and
`ARKStepSetMaxStepsBetweenJac` have been removed and replaced with
`ARKStepSetLSetupFrequency` and `ARKStepSetMaxStepsBetweenJac` respectively.

The previously deprecated function `CVodeSetMaxStepsBetweenJac` has been removed
and replaced with `CVodeSetJacEvalFrequency`.

The ARKODE, CVODE, IDA, and KINSOL Fortran 77 interfaces have been removed. See
the "SUNDIALS Fortran Interface" section in the user guides and the F2003
example programs for more details using the SUNDIALS Fortran 2003 module
interfaces.

#### Namespace Changes

The CUDA, HIP, and SYCL execution policies have been moved from the `sundials`
namespace to the `sundials::cuda`, `sundials::hip`, and `sundials::sycl`
namespaces respectively. Accordingly, the prefixes "Cuda", "Hip", and "Sycl"
have been removed from the execution policy classes and methods.

The `Sundials` namespace used by the Trilinos Tpetra NVector has been replaced
with the `sundials::trilinos::nvector_tpetra` namespace.

### Major Features

#### SUNProfiler

A capability to profile/instrument SUNDIALS library code has been added. This
can be enabled with the CMake option `SUNDIALS_BUILD_WITH_PROFILING`. A built-in
profiler will be used by default, but the
[Caliper](https://github.com/LLNL/Caliper) library can also be used instead with
the CMake option `ENABLE_CALIPER`. See the documentation section on profiling
for more details. **WARNING**: Profiling will impact performance, and should be
enabled judiciously.

#### IMEX MRI Methods and MRIStepInnerStepper Object

The ARKODE MRIStep module has been extended to support implicit-explicit (IMEX)
multirate infinitesimal generalized additive Runge-Kutta (MRI-GARK) methods. As
such, `MRIStepCreate` has been updated to include arguments for the slow
explicit and slow implicit ODE right-hand side functions. `MRIStepCreate` has
also been updated to require attaching an `MRIStepInnerStepper` for evolving the
fast time scale. `MRIStepReInit` has been similarly updated to take explicit
and implicit right-hand side functions as input. Codes using explicit or
implicit MRI methods will need to update `MRIStepCreate` and `MRIStepReInit`
calls to pass `NULL` for either the explicit or implicit right-hand side
function as appropriate. If ARKStep is used as the fast time scale integrator,
codes will need to call `ARKStepCreateMRIStepInnerStepper` to wrap the ARKStep
memory as an `MRIStepInnerStepper` object. Additionally, `MRIStepGetNumRhsEvals`
has been updated to return the number of slow implicit and explicit function
evaluations. The coupling table structure `MRIStepCouplingMem` and the
functions `MRIStepCoupling_Alloc` and `MRIStepCoupling_Create` have also
been updated to support IMEX-MRI-GARK methods.

### New Features

Two new optional vector operations, `N_VDotProdMultiLocal` and
`N_VDotProdMultiAllReduce`, have been added to support low-synchronization
methods for Anderson acceleration.

The implementation of solve-decoupled implicit MRI-GARK methods has been updated
to remove extraneous slow implicit function calls and reduce the memory
requirements.

Added a new function `CVodeGetLinSolveStats` to get the CVODES linear solver
statistics as a group.

Added a new function, `CVodeSetMonitorFn`, that takes a user-function
to be called by CVODES after every `nst` successfully completed time-steps.
This is intended to provide a way of monitoring the CVODES statistics
throughout the simulation.

New orthogonalization methods were added for use within Anderson acceleration
in KINSOL. See the "Anderson Acceleration QR Factorization" subsection within
the mathematical considerations chapter of the user guide and the
`KINSetOrthAA` function documentation for more details.

### Deprecation Notice

The serial, PThreads, PETSc, *hypre*, Parallel, OpenMP_DEV, and OpenMP vector
functions `N_VCloneVectorArray_*` and `N_VDestroyVectorArray_*` have been
deprecated. The generic `N_VCloneVectorArray` and `N_VDestroyVectorArray`
functions should be used instead.

Many constants, types, and functions have been renamed so that they are properly
namespaced. The old names have been deprecated and will be removed in SUNDIALS
v7.0.0.

The following constants, macros, and typedefs are now deprecated:

| Deprecated Name            | New Name                          |
|:---------------------------|:----------------------------------|
| `realtype`                 | `sunrealtype`                     |
| `booleantype`              | `sunbooleantype`                  |
| `RCONST`                   | `SUN_RCONST`                      |
| `BIG_REAL`                 | `SUN_BIG_REAL`                    |
| `SMALL_REAL`               | `SUN_SMALL_REAL`                  |
| `UNIT_ROUNDOFF`            | `SUN_UNIT_ROUNDOFF`               |
| `PREC_NONE`                | `SUN_PREC_NONE`                   |
| `PREC_LEFT`                | `SUN_PREC_LEFT`                   |
| `PREC_RIGHT`               | `SUN_PREC_RIGHT`                  |
| `PREC_BOTH`                | `SUN_PREC_BOTH`                   |
| `MODIFIED_GS`              | `SUN_MODIFIED_GS`                 |
| `CLASSICAL_GS`             | `SUN_CLASSICAL_GS`                |
| `ATimesFn`                 | `SUNATimesFn`                     |
| `PSetupFn`                 | `SUNPSetupFn`                     |
| `PSolveFn`                 | `SUNPSolveFn`                     |
| `DlsMat`                   | `SUNDlsMat`                       |
| `DENSE_COL`                | `SUNDLS_DENSE_COL`                |
| `DENSE_ELEM`               | `SUNDLS_DENSE_ELEM`               |
| `BAND_COL`                 | `SUNDLS_BAND_COL`                 |
| `BAND_COL_ELEM`            | `SUNDLS_BAND_COL_ELEM`            |
| `BAND_ELEM`                | `SUNDLS_BAND_ELEM`                |
| `SDIRK_2_1_2`              | `ARKODE_SDIRK_2_1_2`              |
| `BILLINGTON_3_3_2`         | `ARKODE_BILLINGTON_3_3_2`         |
| `TRBDF2_3_3_2`             | `ARKODE_TRBDF2_3_3_2`             |
| `KVAERNO_4_2_3`            | `ARKODE_KVAERNO_4_2_3`            |
| `ARK324L2SA_DIRK_4_2_3`    | `ARKODE_ARK324L2SA_DIRK_4_2_3`    |
| `CASH_5_2_4`               | `ARKODE_CASH_5_2_4`               |
| `CASH_5_3_4`               | `ARKODE_CASH_5_3_4`               |
| `SDIRK_5_3_4`              | `ARKODE_SDIRK_5_3_4`              |
| `KVAERNO_5_3_4`            | `ARKODE_KVAERNO_5_3_4`            |
| `ARK436L2SA_DIRK_6_3_4`    | `ARKODE_ARK436L2SA_DIRK_6_3_4`    |
| `KVAERNO_7_4_5`            | `ARKODE_KVAERNO_7_4_5`            |
| `ARK548L2SA_DIRK_8_4_5`    | `ARKODE_ARK548L2SA_DIRK_8_4_5`    |
| `ARK437L2SA_DIRK_7_3_4`    | `ARKODE_ARK437L2SA_DIRK_7_3_4`    |
| `ARK548L2SAb_DIRK_8_4_5`   | `ARKODE_ARK548L2SAb_DIRK_8_4_5`   |
| `MIN_DIRK_NUM`             | `ARKODE_MIN_DIRK_NUM`             |
| `MAX_DIRK_NUM`             | `ARKODE_MAX_DIRK_NUM`             |
| `MIS_KW3`                  | `ARKODE_MIS_KW3`                  |
| `MRI_GARK_ERK33a`          | `ARKODE_MRI_GARK_ERK33a`          |
| `MRI_GARK_ERK45a`          | `ARKODE_MRI_GARK_ERK45a`          |
| `MRI_GARK_IRK21a`          | `ARKODE_MRI_GARK_IRK21a`          |
| `MRI_GARK_ESDIRK34a`       | `ARKODE_MRI_GARK_ESDIRK34a`       |
| `MRI_GARK_ESDIRK46a`       | `ARKODE_MRI_GARK_ESDIRK46a`       |
| `IMEX_MRI_GARK3a`          | `ARKODE_IMEX_MRI_GARK3a`          |
| `IMEX_MRI_GARK3b`          | `ARKODE_IMEX_MRI_GARK3b`          |
| `IMEX_MRI_GARK4`           | `ARKODE_IMEX_MRI_GARK4`           |
| `MIN_MRI_NUM`              | `ARKODE_MIN_MRI_NUM`              |
| `MAX_MRI_NUM`              | `ARKODE_MAX_MRI_NUM`              |
| `DEFAULT_MRI_TABLE_3`      | `MRISTEP_DEFAULT_TABLE_3`         |
| `DEFAULT_EXPL_MRI_TABLE_3` | `MRISTEP_DEFAULT_EXPL_TABLE_3`    |
| `DEFAULT_EXPL_MRI_TABLE_4` | `MRISTEP_DEFAULT_EXPL_TABLE_4`    |
| `DEFAULT_IMPL_SD_TABLE_2`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_2` |
| `DEFAULT_IMPL_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_3` |
| `DEFAULT_IMPL_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMPL_SD_TABLE_4` |
| `DEFAULT_IMEX_SD_TABLE_3`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_3` |
| `DEFAULT_IMEX_SD_TABLE_4`  | `MRISTEP_DEFAULT_IMEX_SD_TABLE_4` |
| `HEUN_EULER_2_1_2`         | `ARKODE_HEUN_EULER_2_1_2`         |
| `BOGACKI_SHAMPINE_4_2_3`   | `ARKODE_BOGACKI_SHAMPINE_4_2_3`   |
| `ARK324L2SA_ERK_4_2_3`     | `ARKODE_ARK324L2SA_ERK_4_2_3`     |
| `ZONNEVELD_5_3_4`          | `ARKODE_ZONNEVELD_5_3_4`          |
| `ARK436L2SA_ERK_6_3_4`     | `ARKODE_ARK436L2SA_ERK_6_3_4`     |
| `SAYFY_ABURUB_6_3_4`       | `ARKODE_SAYFY_ABURUB_6_3_4`       |
| `CASH_KARP_6_4_5`          | `ARKODE_CASH_KARP_6_4_5`          |
| `FEHLBERG_6_4_5`           | `ARKODE_FEHLBERG_6_4_5`           |
| `DORMAND_PRINCE_7_4_5`     | `ARKODE_DORMAND_PRINCE_7_4_5`     |
| `ARK548L2SA_ERK_8_4_5`     | `ARKODE_ARK548L2SA_ERK_8_4_5`     |
| `VERNER_8_5_6`             | `ARKODE_VERNER_8_5_6`             |
| `FEHLBERG_13_7_8`          | `ARKODE_FEHLBERG_13_7_8`          |
| `KNOTH_WOLKE_3_3`          | `ARKODE_KNOTH_WOLKE_3_3`          |
| `ARK437L2SA_ERK_7_3_4`     | `ARKODE_ARK437L2SA_ERK_7_3_4`     |
| `ARK548L2SAb_ERK_8_4_5`    | `ARKODE_ARK548L2SAb_ERK_8_4_5`    |
| `MIN_ERK_NUM`              | `ARKODE_MIN_ERK_NUM`              |
| `MAX_ERK_NUM`              | `ARKODE_MAX_ERK_NUM`              |
| `DEFAULT_ERK_2`            | `ARKSTEP_DEFAULT_ERK_2`           |
| `DEFAULT_ERK_3`            | `ARKSTEP_DEFAULT_ERK_3`           |
| `DEFAULT_ERK_4`            | `ARKSTEP_DEFAULT_ERK_4`           |
| `DEFAULT_ERK_5`            | `ARKSTEP_DEFAULT_ERK_5`           |
| `DEFAULT_ERK_6`            | `ARKSTEP_DEFAULT_ERK_6`           |
| `DEFAULT_ERK_8`            | `ARKSTEP_DEFAULT_ERK_8`           |
| `DEFAULT_DIRK_2`           | `ARKSTEP_DEFAULT_DIRK_2`          |
| `DEFAULT_DIRK_3`           | `ARKSTEP_DEFAULT_DIRK_3`          |
| `DEFAULT_DIRK_4`           | `ARKSTEP_DEFAULT_DIRK_4`          |
| `DEFAULT_DIRK_5`           | `ARKSTEP_DEFAULT_DIRK_5`          |
| `DEFAULT_ARK_ETABLE_3`     | `ARKSTEP_DEFAULT_ARK_ETABLE_3`    |
| `DEFAULT_ARK_ETABLE_4`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ETABLE_5`     | `ARKSTEP_DEFAULT_ARK_ETABLE_4`    |
| `DEFAULT_ARK_ITABLE_3`     | `ARKSTEP_DEFAULT_ARK_ITABLE_3`    |
| `DEFAULT_ARK_ITABLE_4`     | `ARKSTEP_DEFAULT_ARK_ITABLE_4`    |
| `DEFAULT_ARK_ITABLE_5`     | `ARKSTEP_DEFAULT_ARK_ITABLE_5`    |
| `DEFAULT_ERK_2`            | `ERKSTEP_DEFAULT_2`               |
| `DEFAULT_ERK_3`            | `ERKSTEP_DEFAULT_3`               |
| `DEFAULT_ERK_4`            | `ERKSTEP_DEFAULT_4`               |
| `DEFAULT_ERK_5`            | `ERKSTEP_DEFAULT_5`               |
| `DEFAULT_ERK_6`            | `ERKSTEP_DEFAULT_6`               |
| `DEFAULT_ERK_8`            | `ERKSTEP_DEFAULT_8`               |

In addition, the following functions are now deprecated (compile-time warnings
will be printed if supported by the compiler):

| Deprecated Name               | New Name                     |
|:------------------------------|:-----------------------------|
| `CVSpilsSetLinearSolver`      | `CVodeSetLinearSolver`       |
| `CVSpilsSetEpsLin`            | `CVodeSetEpsLin`             |
| `CVSpilsSetPreconditioner`    | `CVodeSetPreconditioner`     |
| `CVSpilsSetJacTimes`          | `CVodeSetJacTimes`           |
| `CVSpilsGetWorkSpace`         | `CVodeGetLinWorkSpace`       |
| `CVSpilsGetNumPrecEvals`      | `CVodeGetNumPrecEvals`       |
| `CVSpilsGetNumPrecSolves`     | `CVodeGetNumPrecSolves`      |
| `CVSpilsGetNumLinIters`       | `CVodeGetNumLinIters`        |
| `CVSpilsGetNumConvFails`      | `CVodeGetNumConvFails`       |
| `CVSpilsGetNumJTSetupEvals`   | `CVodeGetNumJTSetupEvals`    |
| `CVSpilsGetNumJtimesEvals`    | `CVodeGetNumJtimesEvals`     |
| `CVSpilsGetNumRhsEvals`       | `CVodeGetNumLinRhsEvals`     |
| `CVSpilsGetLastFlag`          | `CVodeGetLastLinFlag`        |
| `CVSpilsGetReturnFlagName`    | `CVodeGetLinReturnFlagName`  |
| `CVSpilsSetLinearSolverB`     | `CVodeSetLinearSolverB`      |
| `CVSpilsSetEpsLinB`           | `CVodeSetEpsLinB`            |
| `CVSpilsSetPreconditionerB`   | `CVodeSetPreconditionerB`    |
| `CVSpilsSetPreconditionerBS`  | `CVodeSetPreconditionerBS`   |
| `CVSpilsSetJacTimesB`         | `CVodeSetJacTimesB`          |
| `CVSpilsSetJacTimesBS`        | `CVodeSetJacTimesBS`         |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `CVDlsSetLinearSolverB`       | `CVodeSetLinearSolverB`      |
| `CVDlsSetJacFnB`              | `CVodeSetJacFnB`             |
| `CVDlsSetJacFnBS`             | `CVodeSetJacFnBS`            |
| `CVDlsSetLinearSolver`        | `CVodeSetLinearSolver`       |
| `CVDlsSetJacFn`               | `CVodeSetJacFn`              |
| `CVDlsGetWorkSpace`           | `CVodeGetLinWorkSpace`       |
| `CVDlsGetNumJacEvals`         | `CVodeGetNumJacEvals`        |
| `CVDlsGetNumRhsEvals`         | `CVodeGetNumLinRhsEvals`     |
| `CVDlsGetLastFlag`            | `CVodeGetLastLinFlag`        |
| `CVDlsGetReturnFlagName`      | `CVodeGetLinReturnFlagName`  |
| `KINDlsSetLinearSolver`       | `KINSetLinearSolver`         |
| `KINDlsSetJacFn`              | `KINSetJacFn`                |
| `KINDlsGetWorkSpace`          | `KINGetLinWorkSpace`         |
| `KINDlsGetNumJacEvals`        | `KINGetNumJacEvals`          |
| `KINDlsGetNumFuncEvals`       | `KINGetNumLinFuncEvals`      |
| `KINDlsGetLastFlag`           | `KINGetLastLinFlag`          |
| `KINDlsGetReturnFlagName`     | `KINGetLinReturnFlagName`    |
| `KINSpilsSetLinearSolver`     | `KINSetLinearSolver`         |
| `KINSpilsSetPreconditioner`   | `KINSetPreconditioner`       |
| `KINSpilsSetJacTimesVecFn`    | `KINSetJacTimesVecFn`        |
| `KINSpilsGetWorkSpace`        | `KINGetLinWorkSpace`         |
| `KINSpilsGetNumPrecEvals`     | `KINGetNumPrecEvals`         |
| `KINSpilsGetNumPrecSolves`    | `KINGetNumPrecSolves`        |
| `KINSpilsGetNumLinIters`      | `KINGetNumLinIters`          |
| `KINSpilsGetNumConvFails`     | `KINGetNumLinConvFails`      |
| `KINSpilsGetNumJtimesEvals`   | `KINGetNumJtimesEvals`       |
| `KINSpilsGetNumFuncEvals`     | `KINGetNumLinFuncEvals`      |
| `KINSpilsGetLastFlag`         | `KINGetLastLinFlag`          |
| `KINSpilsGetReturnFlagName`   | `KINGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolver`     | `IDASetLinearSolver`         |
| `IDASpilsSetPreconditioner`   | `IDASetPreconditioner`       |
| `IDASpilsSetJacTimes`         | `IDASetJacTimes`             |
| `IDASpilsSetEpsLin`           | `IDASetEpsLin`               |
| `IDASpilsSetIncrementFactor`  | `IDASetIncrementFactor`      |
| `IDASpilsGetWorkSpace`        | `IDAGetLinWorkSpace`         |
| `IDASpilsGetNumPrecEvals`     | `IDAGetNumPrecEvals`         |
| `IDASpilsGetNumPrecSolves`    | `IDAGetNumPrecSolves`        |
| `IDASpilsGetNumLinIters`      | `IDAGetNumLinIters`          |
| `IDASpilsGetNumConvFails`     | `IDAGetNumLinConvFails`      |
| `IDASpilsGetNumJTSetupEvals`  | `IDAGetNumJTSetupEvals`      |
| `IDASpilsGetNumJtimesEvals`   | `IDAGetNumJtimesEvals`       |
| `IDASpilsGetNumResEvals`      | `IDAGetNumLinResEvals`       |
| `IDASpilsGetLastFlag`         | `IDAGetLastLinFlag`          |
| `IDASpilsGetReturnFlagName`   | `IDAGetLinReturnFlagName`    |
| `IDASpilsSetLinearSolverB`    | `IDASetLinearSolverB`        |
| `IDASpilsSetEpsLinB`          | `IDASetEpsLinB`              |
| `IDASpilsSetIncrementFactorB` | `IDASetIncrementFactorB`     |
| `IDASpilsSetPreconditionerB`  | `IDASetPreconditionerB`      |
| `IDASpilsSetPreconditionerBS` | `IDASetPreconditionerBS`     |
| `IDASpilsSetJacTimesB`        | `IDASetJacTimesB`            |
| `IDASpilsSetJacTimesBS`       | `IDASetJacTimesBS`           |
| `IDADlsSetLinearSolver`       | `IDASetLinearSolver`         |
| `IDADlsSetJacFn`              | `IDASetJacFn`                |
| `IDADlsGetWorkSpace`          | `IDAGetLinWorkSpace`         |
| `IDADlsGetNumJacEvals`        | `IDAGetNumJacEvals`          |
| `IDADlsGetNumResEvals`        | `IDAGetNumLinResEvals`       |
| `IDADlsGetLastFlag`           | `IDAGetLastLinFlag`          |
| `IDADlsGetReturnFlagName`     | `IDAGetLinReturnFlagName`    |
| `IDADlsSetLinearSolverB`      | `IDASetLinearSolverB`        |
| `IDADlsSetJacFnB`             | `IDASetJacFnB`               |
| `IDADlsSetJacFnBS`            | `IDASetJacFnBS`              |
| `DenseGETRF`                  | `SUNDlsMat_DenseGETRF`       |
| `DenseGETRS`                  | `SUNDlsMat_DenseGETRS`       |
| `denseGETRF`                  | `SUNDlsMat_denseGETRF`       |
| `denseGETRS`                  | `SUNDlsMat_denseGETRS`       |
| `DensePOTRF`                  | `SUNDlsMat_DensePOTRF`       |
| `DensePOTRS`                  | `SUNDlsMat_DensePOTRS`       |
| `densePOTRF`                  | `SUNDlsMat_densePOTRF`       |
| `densePOTRS`                  | `SUNDlsMat_densePOTRS`       |
| `DenseGEQRF`                  | `SUNDlsMat_DenseGEQRF`       |
| `DenseORMQR`                  | `SUNDlsMat_DenseORMQR`       |
| `denseGEQRF`                  | `SUNDlsMat_denseGEQRF`       |
| `denseORMQR`                  | `SUNDlsMat_denseORMQR`       |
| `DenseCopy`                   | `SUNDlsMat_DenseCopy`        |
| `denseCopy`                   | `SUNDlsMat_denseCopy`        |
| `DenseScale`                  | `SUNDlsMat_DenseScale`       |
| `denseScale`                  | `SUNDlsMat_denseScale`       |
| `denseAddIdentity`            | `SUNDlsMat_denseAddIdentity` |
| `DenseMatvec`                 | `SUNDlsMat_DenseMatvec`      |
| `denseMatvec`                 | `SUNDlsMat_denseMatvec`      |
| `BandGBTRF`                   | `SUNDlsMat_BandGBTRF`        |
| `bandGBTRF`                   | `SUNDlsMat_bandGBTRF`        |
| `BandGBTRS`                   | `SUNDlsMat_BandGBTRS`        |
| `bandGBTRS`                   | `SUNDlsMat_bandGBTRS`        |
| `BandCopy`                    | `SUNDlsMat_BandCopy`         |
| `bandCopy`                    | `SUNDlsMat_bandCopy`         |
| `BandScale`                   | `SUNDlsMat_BandScale`        |
| `bandScale`                   | `SUNDlsMat_bandScale`        |
| `bandAddIdentity`             | `SUNDlsMat_bandAddIdentity`  |
| `BandMatvec`                  | `SUNDlsMat_BandMatvec`       |
| `bandMatvec`                  | `SUNDlsMat_bandMatvec`       |
| `ModifiedGS`                  | `SUNModifiedGS`              |
| `ClassicalGS`                 | `SUNClassicalGS`             |
| `QRfact`                      | `SUNQRFact`                  |
| `QRsol`                       | `SUNQRsol`                   |
| `DlsMat_NewDenseMat`          | `SUNDlsMat_NewDenseMat`      |
| `DlsMat_NewBandMat`           | `SUNDlsMat_NewBandMat`       |
| `DestroyMat`                  | `SUNDlsMat_DestroyMat`       |
| `NewIntArray`                 | `SUNDlsMat_NewIntArray`      |
| `NewIndexArray`               | `SUNDlsMat_NewIndexArray`    |
| `NewRealArray`                | `SUNDlsMat_NewRealArray`     |
| `DestroyArray`                | `SUNDlsMat_DestroyArray`     |
| `AddIdentity`                 | `SUNDlsMat_AddIdentity`      |
| `SetToZero`                   | `SUNDlsMat_SetToZero`        |
| `PrintMat`                    | `SUNDlsMat_PrintMat`         |
| `newDenseMat`                 | `SUNDlsMat_newDenseMat`      |
| `newBandMat`                  | `SUNDlsMat_newBandMat`       |
| `destroyMat`                  | `SUNDlsMat_destroyMat`       |
| `newIntArray`                 | `SUNDlsMat_newIntArray`      |
| `newIndexArray`               | `SUNDlsMat_newIndexArray`    |
| `newRealArray`                | `SUNDlsMat_newRealArray`     |
| `destroyArray`                | `SUNDlsMat_destroyArray`     |

In addition, the entire `sundials_lapack.h` header file is now deprecated for
removal in SUNDIALS v7.0.0. Note, this header file is not needed to use the
SUNDIALS LAPACK linear solvers.

Deprecated ARKODE nonlinear solver predictors: specification of the ARKStep
"bootstrap" or "minimum correction" predictors (options 4 and 5 from
`ARKStepSetPredictorMethod`), or MRIStep "bootstrap" predictor (option 4 from
`MRIStepSetPredictorMethod`), will output a deprecation warning message.
These options will be removed in a future release.

