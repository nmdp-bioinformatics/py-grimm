import numpy as np
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize


with open("README.md") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read().split("\n")

setup(
    name="GRMA",
    version="0.0.1",
    author="Amit Kabya",
    author_email="kabya.amit@gmail.com",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    description="Graph Based Matching",
    install_requires=requirements,
    license="LGPL 3.0",
    long_description=readme,
    long_description_content_type="text/markdown",
    include_package_data=True,
    keywords="grma",
    scripts=["scripts/build-imputation-validation.sh", "scripts/runfile.py"],
    packages=find_packages(
        include=[
            "GRMA",
            "GRMA.Build",
            "GRMA.Grim",
            "GRMA.Match",
            "GRMA.Utilities"
        ]
    ),
    test_suite="tests",
    zip_safe=False,
    ext_modules=cythonize(
        [
            Extension("GRMA.Utilities.cutils", ["GRMA/Utilities/cutils.pyx"],
                      define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]),
            Extension("GRMA.Match.lol_graph", ["GRMA/Match/lol_graph.pyx"],
                      define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")])
        ]
    ),
    include_dirs=[np.get_include()]
)
