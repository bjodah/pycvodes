from skbuild import setup  # This line replaces 'from setuptools import setup'
from setuptools import find_packages

setup(
    name="pycvodes",
    version="0.14.7",
    description="A python wrapper around SUNDIALS' CVODES",
    author='Björn Dahlgren',
    license="BSD-2-Clause",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=["numpy"],
)
