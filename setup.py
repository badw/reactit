
"""
reactit: Reaction Iterator
"""

from os.path import abspath, dirname
from setuptools import find_packages, setup

setup(
    name='reactit',
    version='0.0.1',
    description='reaction iterator',
    url="https://github.com/badw/reactit",
    author="Benjamin A. D. Williamson",
    author_email="benjamin.williamson@ntnu.no",
    license='MIT',
    packages=['reactit'],
    install_requires=[
        'numpy',
        'chempy',
        'tqdm',
        'pathos',
        ],
    python_requires=">=3.11",
    classifiers=[
            "Programming Language :: Python :: 3.11",
            "Development Status :: 3 - Beta",
            "Intended Audience :: Science/Research",
            "Intended Audience :: System Administrators",
            "Intended Audience :: Industry",
            "Intended Audience :: Information Technology",
            "Operating System :: OS Independent",
            "Topic :: Other/Nonlisted Topic",
            "Topic :: Scientific/Engineering",
        ],
    #entry_points = {
    #    "console_scripts":[
    #        "arcs-app = app.arcs_app:start"
    #        ]
    #    },
    #include_package_data=True
    )
