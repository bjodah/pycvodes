cmake_minimum_required(VERSION 3.15)
project(pycvodes LANGUAGES CXX)

# Find Python and SUNDIALS
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(SUNDIALS REQUIRED)
find_package(NumPy REQUIRED)
find_package(Cython REQUIRED)
include(UseCython)

# Create the extension
cython_transpile(pycvodes/_cvodes.pyx LANGUAGE CXX OUTPUT_VARIABLE _cvodes_cxx)
python_add_library(_cvodes MODULE "${_cvodes_cxx}" WITH_SOABI)

target_link_libraries(_cvodes PRIVATE
    ${SUNDIALS_LIBRARIES}
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
