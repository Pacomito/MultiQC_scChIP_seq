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
log.info("IN CUSTOM_CODE.PY")
# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.scChIPseq_version = get_distribution("scChIPseq").version


# Add default config options for the things that are used in MultiQC_NGI
def scChIPseq_plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.
    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None

    log.info("Running scChIPseq MultiQC Plugin v{}".format(config.scChIPseq_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    # Add to the search patterns used by modules
    if 'scChIPseq/all_logs' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/all_logs': { 'fn': '*scChIPseq_logs.txt' } } )
        log.info("Added scChIPseq/all_logs to the search patterns")
    if 'scChIPseq/flagged_count' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/flagged_count': { 'fn': '*_flagged.count' } } )
        log.info("Added scChIPseq/flagged_count to the search patterns")
    if 'scChIPseq/flagged_PCR_count' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/flagged_PCR_count': { 'fn': '*_flagged_rmPCR.count' } } )
        log.info("Added scChIPseq/flagged_PCR_count to the search patterns")
    if 'scChIPseq/flagged_PCR_RT_count' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/flagged_PCR_RT_count': { 'fn': '*_flagged_rmPCR_RT.count' } } )
        log.info("Added scChIPseq/flagged_PCR_RT_count to the search patterns")
    if 'scChIPseq/flagged_PCR_RT_rmDup_count' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/flagged_PCR_RT_rmDup_count': { 'fn': '*_rmDup.count' } } )
        log.info("Added scChIPseq/flagged_PCR_RT_rmDup_count to the search patterns")
    if 'scChIPseq/count_matrix' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/count_matrix': { 'fn': '*.tsv' } } )
        log.info("Added scChIPseq/count_matrix to the search patterns")
    if 'scChIPseq/bigwig:' not in config.sp:
        config.update_dict( config.sp, { 'scChIPseq/bigwig:': { 'fn': '*_rmDup.bw' } } )
        log.info("Added scChIPseq/bigwig to the search patterns")
