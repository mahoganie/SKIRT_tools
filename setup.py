import setuptools
from src.skirt_tools import __version__

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SKIRT_tools",
    version=__version__,
    author="Jenny Wang",
    author_email="jennydw@uci.edu",
    description="Tool to run SKIRT simulations and make mock images.",
    url="https://github.com/mahoganie/SKIRT-tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'astropy',
        'firestudio',
        'scipy',
        'pandas',
        'h5py',
        "lxml>=5.3.0",
        "pathlib>=1.0.1",
        "numpy<=2.0.0",
        "pandas<=1.1.5",
        "matplotlib<=3.3.4",
        "abg_python",
    ],
)
