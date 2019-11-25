# -*- coding: utf-8 -*-

import sys
import os
import logging
from glob import iglob
from collections import Counter, OrderedDict

def readable_file(fn):
    """Check if the file is readable"""
    fn = os.path.abspath(fn)
    if not os.path.isfile (fn) or not os.access (fn, os.R_OK):
        raise ont2cramError("File '{}' does not exist or is not readable".format(fn))

def writable_dir(fn):
    """Check if the dir is writable"""
    fn = os.path.abspath(fn)
    if not os.path.isdir(fn) or not os.access (fn, os.W_OK):
        raise ont2cramError("Directory '{}' does not exist or is not writable".format(fn))

def readable_dir(fn):
    """Check if the dir is readable"""
    fn = os.path.abspath(fn)
    if not os.path.isdir(fn) or not os.access (fn, os.R_OK):
        raise ont2cramError("Directory '{}' does not exist or is not readable".format(fn))

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise ont2cramError ("Error creating output folder '{}'".format(fn))

def get_logger (name=None, verbose=False, quiet=False):
    """Set logger to appropriate log level"""

    logging.basicConfig(format='%(message)s')
    logger = logging.getLogger(name)
    # Define overall verbose level
    if verbose:
        logger.setLevel(logging.DEBUG)
    elif quiet:
        logger.setLevel(logging.WARNING)
    else:
        logger.setLevel(logging.INFO)
    return logger

def recursive_file_gen (dir, ext):
    """
    create a generator listing all files with a particular extension in a folder arborescence
    The recursivity is broken when at least 1 file with a particular extenssion is found.
    """
    # In the case where the folder is a file
    if os.path.isdir(dir):

        # If matching files in the folder
        file_found=False
        for fn in iglob (os.path.join(dir, "*."+ext)):
            yield fn
            file_found=True

        # If no matching file go deeper until a leaf containing fast5 is found
        if not file_found:
            for item in listdir(dir):
                for fn in recursive_file_gen (os.path.join(dir, item), ext):
                    yield fn

def dict_to_str (d, tab="\t", ntab=0):
    """ Transform a multilevel dict to a tabulated str """
    m = ""

    if isinstance(d, Counter):
        for i, j in d.most_common():
            m += "{}{}: {:,}\n".format(tab*ntab, i, j)

    else:
        for i, j in d.items():
            if isinstance(j, dict):
                j = dict_to_str(j, tab=tab, ntab=ntab+1)
                m += "{}{}\n{}".format(tab*ntab, i, j)
            else:
                m += "{}{}: {}\n".format(tab*ntab, i, j)
    return m

class ont2cramError (Exception):
    """ Basic exception class for ont2cram package """
    pass
