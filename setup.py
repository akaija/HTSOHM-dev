#!/usr/bin/env python

from distutils.core import setup
setup(
    name = 'HTSOHM',
    version = '0.1.0',
    description = 'High-throughput Screening of Hypothetical Materials',
    author = 'A Reino Kaija',
    author_email = 'ark111@pitt.edu',
    url = 'https://github.com/akaija/HTSOHM-dev',
    packages = ['htsohm', 'htsohm.db', 'htsohm.simulation',
                'htsohm.pseudomaterial_generator']
)
