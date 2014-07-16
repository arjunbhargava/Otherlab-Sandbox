#!/usr/bin/env python

from setuptools import setup,find_packages

setup(
  # Basics
  name='ompl_planning',
  version='0.0-dev',
  description='ompl path planning',
  author='Otherlab et al.',
  author_email='arjunb@stanford.edu',
  url='',

  # Installation
  packages=find_packages(),
  package_data={'kinematic_chain':['*.py','*.so']},
)
