from setuptools import setup
import setuptools


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='BiofilmSimulation',
    version='1.0.0',
    description='The model applies methods from molecular dynamics (MD) and takes into account different '
                'physical and biological effects. '
                'The software provides great flexibility by enabling the user to switch easily between sets of '
                'constants e.g. to model different bacterial strains. '
                'Furthermore, the software includes functions for visualisation of the models behaviour over time.',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/StudDavid/biofilm_growth_modeling',
    author='David Theidel, Thorben Klamt, Tizian Dege, Jonas Müller',
    author_email='theidel@stud.uni-hannover.de',
    license='MIT License',
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)

