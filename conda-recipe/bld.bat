> pycvodes/_config.py (
  @echo.env = {
  @echo.  "LAPACK": "",
  @echo.  "SUNDIALS_LIBS": "sundials_cvodes,sundials_nvecserial,sundials_sunlinsoldense,sundials_sunlinsolband,sundials_sunlinsolspgmr,sundials_sunlinsolspbcgs,sundials_sunlinsolsptfqmr,sundials_sunmatrixdense,sundials_sunmatrixband",
  @echo.  "NO_KLU": "1",
  @echo.  "NO_LAPACK": "1",
  @echo.  "SUNDIALS_PRECISION": "double",
  @echo.  "REAL_TYPE": "double",
  @echo.  "INDEX_TYPE": "int32_t"
  @echo.}
)
python -m pip install --no-deps --ignore-installed .
