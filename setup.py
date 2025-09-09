from skbuild import setup  # This line replaces 'from setuptools import setup'
setup(
    name="pycvodes",
    version="0.14.7",
    description="A python wrapper around SUNDIALS' CVODES",
    author='Björn Dahlgren',
    license="BSD-2-Clause",
    packages=['pycvodes'],
    python_requires=">=3.8",
    install_requires=["numpy"],
)
