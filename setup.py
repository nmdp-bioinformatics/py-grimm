import numpy as np
from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize


with open("grma/README.md") as readme_file:
    readme = readme_file.read()

with open("requirements.txt") as requirements_file:
    requirements = requirements_file.read().split("\n")

setup(
    name="GRMA",
    version="0.0.1",
    author="Amit Kabya",
    maintainer="Ziv Naim",
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
    keywords="graph,em,family",
    packages=find_packages(
        include=[
            "gram",
            "grem",
            "grem.conf",
            "grem.data",
            "grem.EM",
            "grma",
            "grma.donorsgraph",
            "grma.imputation",
            "grma.match",
            "grma.utilities"
        ]
    ),
    test_suite="tests",
    zip_safe=False,
    ext_modules=cythonize(
        [
            Extension("grma.utilities.cutils", ["grma/utilities/cutils.pyx"],
                      define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")]),
            Extension("grma.match.lol_graph", ["grma/match/lol_graph.pyx"],
                      define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")])
        ]
    ),
    include_dirs=[np.get_include()]
)
