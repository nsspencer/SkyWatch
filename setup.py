from setuptools import find_packages, setup

setup(
    name="skywatch",
    version="0.1.0",
    url="https://github.com/nsspencer/SkyPath",
    author="Nathan Spencer",
    description="Aerospace/astrodynamics analysis library providing high level interfaces for coordinates, attitude, access determination, and look angle calculations.",
    packages=find_packages(),
    install_requires=["astropy", "numpy", "portion", "scipy", "pymap3d"],
)