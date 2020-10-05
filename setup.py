#!/usr/bin/env python

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as reqs:
    requirements = reqs.read().splitlines()

setuptools.setup(
    name="sample", # Replace with your own username
    version="0.0.1",
    author="Valerie Hayot-Sasson",
    author_email="valeriehayot@gmail.com",
    description="Adaptive learning on HPC",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/radical-collaboration/SAMPLE-ICM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=requirements
    python_requires='>=3.6',
)
