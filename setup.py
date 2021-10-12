#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="das",
    version="0.1",
    packages=find_packages(exclude=["bin"]),
    scripts=[
        "bin/das",
        "bin/das_convert",
        "bin/das_vasprun2mtp",
        "bin/das_mean_forces",
        "bin/das_atomic_dataset_info",
        "bin/das_diff_atomic_dataset",
        "bin/das_gen_unfitted_pot",
        "bin/das_ambiguity",
    ],
    include_package_data=True,
    python_requires=">=3.7",
    author="hlyang",
    author_email="hlyang1992@gamil.com",
    long_description=long_description,
    description="Dual Adaptive Sampling for Machine Learning Interatomic Potential",
    long_description_content_type="text/markdown",
)
