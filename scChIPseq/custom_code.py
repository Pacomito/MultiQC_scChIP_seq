#!/usr/bin/env python
""" MultiQC example plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')
# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.scChIPseq_version = get_distribution("scChIPseq").version


# Add default config options for the things that are used in MultiQC_NGI
def scChIPseq_plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    log.info("Running scChIPseq MultiQC Plugin v{}".format(config.scChIPseq_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    search_patterns = {
        'scChIPseq/all_logs': { 'fn': '*scChIPseq_logs.txt' },
        'scChIPseq/flagged_count': { 'fn': '*_flagged.count' },
        'scChIPseq/flagged_PCR_count': { 'fn': '*_flagged_rmPCR.count' },
        'scChIPseq/flagged_PCR_RT_count': { 'fn': '*_flagged_rmPCR_RT.count' },
        'scChIPseq/flagged_PCR_RT_rmDup_count': { 'fn': '*_rmDup.count' },
        'scChIPseq/count_matrix': { 'fn': '*.tsv' },
        'scChIPseq/bigwig:': { 'fn': '*_rmDup.bw' }
    }

    # Add to the search patterns used by modules
    for pattern_name, pattern in search_patterns.items():
        if pattern_name not in config.sp:
            config.update_dict( config.sp, { pattern_name: pattern } )
            log.debug("Added {} to the search patterns".format(pattern_name))
