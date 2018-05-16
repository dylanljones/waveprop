# -*- coding: utf-8 -*-
"""
Created on 13 May 2018
@author: Dylan Jones

"""
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
    name='waveprop',
    description='Package for calculating transmission properties of molecular kronig-penney-systems in one dimension',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Dylan Jones',
    author_email='dylanljones94@gmail.com',
    version='1.0.0',
    install_requires=install_requirements,
    packages=['waveprop', 'waveprop.model', 'waveprop.utils', 'waveprop.plotting', 'waveprop.calculation'],
    url='',
    license='MIT',
)
