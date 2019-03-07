#!/usr/bin/env python
"""
Example plugin for MultiQC, showing how to structure code
and plugin hooks to work effectively with the main MultiQC code.
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'scChIPseq',
    version = version,
    author = 'Pacome Prompsy',
    author_email = 'pacome.prompsy@curie.fr',
    description = "MultiQC plugin for single-cell ChIP seq",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://gitlab.curie.fr/data-analysis/ChIP-seq_single-cell_LBC/tree/refonte_R2',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'multiqc',
        'pandas',
        'numpy'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'scChIPseq = scChIPseq.modules.scChIPseq:MultiqcModule',
        ],
	    'multiqc.hooks.v1': [
            'execution_start = scChIPseq.custom_code:scChIPseq_plugin_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
