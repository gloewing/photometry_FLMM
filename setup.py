from setuptools import find_packages, setup
import os

setup(
    name="fast_flmm_rpy2",  # Choose a unique name
    version="0.1.0",
    packages=find_packages(),  # Automatically find all packages
    author="Josh Lawrimore",  # Optional
    author_email="josh.lawrimore@nih.gov",  # Optional
    description="FLMM for analyzing Fiber Photometry data",  # Optional
    long_description=(
        open("README.md").read() if os.path.exists("README.md") else ""
    ),
    long_description_content_type="text/markdown",  # Optional
    install_requires=[
        # List your package dependencies here, for example:
        "rpy2",
        "jupyter",
        "pandas",
        "pytest",
    ],
    python_requires=">=3.9",  # Specify minimum Python version
)
