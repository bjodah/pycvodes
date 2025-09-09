cmake_minimum_required(VERSION 3.15)
project(pycvodes LANGUAGES CXX)

# Find Python and SUNDIALS
find_package(Python3 COMPONENTS Interpreter Development)
find_package(SUNDIALS REQUIRED)

# Create the extension
add_library(_pycvodes MODULE
    pycvodes/_pycvodes.pyx
)

target_link_libraries(_pycvodes 
    ${SUNDIALS_LIBRARIES}
    Python3::Python
)

target_include_directories(_pycvodes PRIVATE 
    ${SUNDIALS_INCLUDE_DIRS}
    ${Python3_INCLUDE_DIRS}
)

# Set properties for Python extension
set_target_properties(_pycvodes PROPERTIES
    PREFIX ""
    SUFFIX "${PYTHON_MODULE_EXTENSION}"
)

# Install
install(TARGETS _pycvodes DESTINATION .)
