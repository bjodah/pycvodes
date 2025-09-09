cmake_minimum_required(VERSION 3.15)
project(pycvodes LANGUAGES CXX)

# Find Python and SUNDIALS
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(Cython REQUIRED)
include(UseCython)

# Create the extension
cython_transpile(
    pycvodes/_cvodes.pyx
    LANGUAGE CXX
    CYTHON_ARGS
    --include-dir "${CMAKE_CURRENT_SOURCE_DIR}/pycvodes/include"
    --include-dir "${CMAKE_CURRENT_SOURCE_DIR}/external/anyode/cython_def"
    OUTPUT_VARIABLE _cvodes_cxx
)  # <-- scikit-build-core + cython + cython-cmake
#add_cython_target(_cvodes pycvodes/_cvodes.pyx CXX PY3 OUTPUT_VAR _cvodes_cxx)  # <-- scikit-build *ONLY*
python_add_library(_cvodes MODULE "${_cvodes_cxx}" WITH_SOABI)

target_link_libraries(_cvodes PRIVATE
    ${SUNDIALS_LIBRARIES}
    Python::NumPy
)

target_include_directories(_cvodes PRIVATE
    ${SUNDIALS_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
    ${NumPy_INCLUDE_DIRS}
    "${CMAKE_CURRENT_SOURCE_DIR}/pycvodes/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/external/anyode/include"
)

# Install
install(TARGETS _cvodes DESTINATION pycvodes)
