from setuptools import setup
from codecs import open as codec_open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with codec_open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Get install requirements from pip-freeze requirements.txt file
with open(path.join(here, 'requirements.txt'), "r") as f:
    lines = f.readlines()
install_requirements = [line[:-1] for line in lines]

setup(
    name='wave_propagation',
    version='1.0.0',
    packages=['data', 'waveprop', 'waveprop.core', 'waveprop.model'],
    url='',
    license='MIT',
    author='Dylan Jones',
    install_requires=["matplotlib", "numpy", "scipy"],
    author_email='',
    description=''
)
