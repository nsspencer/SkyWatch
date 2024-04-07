from setuptools import find_packages, setup

setup(
    name="skywatch",
    version="0.1.0",
    url="https://github.com/nsspencer/SkyWatch",
    author="Nathan Spencer",
    description="Aerospace/astrodynamics analysis library providing high level interfaces for coordinate, attitude, access, and look angle calculations.",
    packages=find_packages(),
    install_requires=["astropy", "numpy", "portion", "scipy", "pymap3d"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
