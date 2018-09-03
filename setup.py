"""Module to setup DivEn to any machine
"""
#!/usr/bin/env python
# -*- coding: utf-8; fill-column: 120 -*-
#
# Copyright (C) 2018 Andrei Duchko <andrey.duchko@gmail.com>


from setuptools import setup

extra_compile_args = []
library_dirs = []

long_description = """Methods for calculating molecular vibrational and rovibrational energy levels
                      using perturbation theories, analysing and resummation (divergent) perturbation series

Developed by Andrei Duchko at the Laboratory of Molecular Spectroscopy at V.E. Zuev Institute of Atmospheric Optics,
Tomsk, Russia.

Original author:    Andrei Duchko <andrey.duchko@gmail.com>
Current maintainer: Andrei Duchko <andrey.duchko@gmail.com>
See the distribution files and THANKS for further contributions.
"""


setup(name="DivEn",
      author              = "Andrei Duchko, LMS laboratory",
      author_email        = "andrey.duchko@gmail.com",
      maintainer          = "Andrei Duchko",
      maintainer_email    = "andrey.duchko@gmail.com",
      url                 = "https://github.com/FuffiKFuffiK/DivEn",
      description         = "Molecular energy spectrum calculator using divergent series resummation techniques",
      version             = "1.0.dev0",
      long_description    = long_description,
      license             = "Open Source",
      test_suite          = 'tests',
     )
