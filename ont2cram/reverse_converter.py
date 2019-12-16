#!/usr/bin/env python3

# Std lib imports
import os
import sys
import re
import shlex

# Third party imports
import h5py
import pysam
import tqdm
import numpy as np
import numpy.lib.recfunctions as rfn

# Local imports
from collections import OrderedDict, Counter
from ont2cram.common import *

# Define global variables
DT_STR_VLEN = h5py.special_dtype(vlen=str)
FILENAME_TAG = "X0"
READ_NUM_TAG_SHORT = "X1"
READ_NUM_TAG_LONG = "X2"
RESERVED_TAGS = [READ_NUM_TAG_SHORT, READ_NUM_TAG_LONG, FILENAME_TAG]
STR_HEX_PATTERN = re.compile(r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$")
COUNTER = Counter()

# Main Function
def reverse_converter (
    input_file,
    output_dir,
    verbose=False,
    quiet=False,
    progress=False,
    **kwargs):
    """CRAM to Fast5 reverse conversion tool"""
    try:
        # Define logger with appropriate verbosity
        logger = get_logger (name="ont2cram_reverse_converter", verbose=verbose, quiet=quiet)

        logger.debug("Check input files")
        readable_file(input_file)
        #writable_dir(output_dir)
        check_destination_exists(input_file, output_dir)

        class Attribute:
          path = ''
          type = ''
          value = None
          is_col = False

        attr_dict = {}
        with pysam.AlignmentFile(input_file, "rc", check_sq=False) as samfile:
            logger.debug("Read CRAM header")
            for comment in samfile.header["CO"]:
                if not comment.startswith(("ATR:","COL:")): continue

                parts = shlex.split(comment)
                part0 = parts[0].split(':')

                a=Attribute()
                a.path   = part0[1]
                a.type   = part0[2]
                a.value  = parts[2][3:] if len(parts)==3 else None
                a.is_col = comment.startswith("COL:")

                tag = parts[1][3:]

                attr_dict[tag]=a

            def is_hex_str( obj, expected_len, is_null_term=False ):
                bytes_to_remove_end = 1 if is_null_term else 0
                return isinstance(obj,bytes) and len(obj)==expected_len and (not is_null_term or obj[-1:]==0) \
                       and STR_HEX_PATTERN.match(obj[:expected_len-bytes_to_remove_end].decode())

            def get_path(hdf_path, read_number_long, read_number_short):
            	if read_number_short:
            		hdf_path = hdf_path.replace("Read_YYY",  read_number_short)
            	if read_number_long:
            		hdf_path = hdf_path.replace("read_XXXXX",read_number_long )
            	return hdf_path

            def get_path_from_dummy(hdf_path, read_number_long, read_number_short):
            	p = get_path(hdf_path, read_number_long, read_number_short)
            	return p.replace("/dummy_attr","")

            def write_hdf_attr(hdf5_file, attr_path, attr_value, attr_type):
                #print(f"path={attr_path}, val={attr_value}, type={attr_type}")
                group_name,_,attr_name =  attr_path.rpartition('/')
                if attr_name=="noname": raise
                try:
                    group = hdf5_file[group_name]
                except KeyError:
                    group = hdf5_file.create_group(group_name)

                try:
                    group.attrs.create(attr_name, attr_value, dtype=attr_type)
                except (TypeError, ValueError):
                    group.attrs[attr_name] = convert_type(attr_value,attr_type)

            logger.debug("CRAM header")
            for read in tqdm.tqdm(samfile.fetch(until_eof=True), unit=" Reads", unit_scale=True, disable=not progress):
                COUNTER["Reads written"]+=1
                fast5_filename = read.get_tag(FILENAME_TAG)
                read_number_long=None
                try:
                	read_number_long   = "read_"+str(read.get_tag(READ_NUM_TAG_LONG))
                except KeyError:
                	pass

                read_number_short=None
                try:
                	read_number_short  = "Read_"+str(read.get_tag(READ_NUM_TAG_SHORT))
                except KeyError:
                	pass

                output_file = os.path.join(output_dir,fast5_filename)
                dir = os.path.dirname(output_file)
                if not os.path.exists(dir) and len(dir)>0 :
                    os.makedirs(dir)

                with h5py.File( output_file, "a" ) as f:
                    if read.query_name != "nofastq":
                        fastq_lines = np.string_("\n".join( [read.query_name, read.query_sequence, '+', pysam.array_to_qualitystring(read.query_qualities)+'\n'] ) )
                        f.create_dataset( "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq", data=fastq_lines )

                    DSETS = {}
                    for tag_name,tag_val in read.get_tags():
                        if tag_name in RESERVED_TAGS: continue
                        a = attr_dict[tag_name]
                        if a.is_col:
                            dset_name,_,col_name = get_path(
                            	a.path,read_number_long,read_number_short ).rpartition('/')

                            if dset_name.endswith("Fastq") and col_name=="noname":
                                continue
                            if dset_name not in DSETS:
                                DSETS[dset_name] = []
                            dset = DSETS[dset_name]

                            if col_name=="noname":
                                dset.append(tag_val)
                            else:
                            	dset.append(
                                    np.array(
                                        list(tag_val.split('\x03')) if a.type.startswith(('S','U')) else tag_val,
                                        dtype=[(col_name, a.type)]
                                    )
                                 )
                    for dset_name,columns in DSETS.items():
                        d = columns[0] if len(columns)==1 else rfn.merge_arrays(columns, flatten=True, usemask=False)
                        f.create_dataset( dset_name, data=d )

                    # write constant values stored in cram header
                    for a in attr_dict.values():
                        if a.is_col:
                            continue
                        if a.value :
                        	write_hdf_attr( f, get_path(a.path,read_number_long, read_number_short), a.value, a.type )

                    # write tags stored in cram records
                    for tag_name, tag_val in read.get_tags():
                        if tag_name in RESERVED_TAGS: continue
                        a = attr_dict[tag_name]
                        if is_empty_hdf_node(a.path):
                        	f.create_group( get_path_from_dummy(a.path,read_number_long,read_number_short) )
                        	continue
                        if a.is_col: continue
                        if a.value != tag_val:
                            write_hdf_attr( f, get_path(a.path,read_number_long, read_number_short), tag_val, a.type )

    finally:
        logger.info(dict_to_str(COUNTER))

# Helper classes and functions
def convert_type(val, typ):
    #print(f"val={val}, typ={typ}, type={type(val)}")
    if typ.startswith('U'):    return str.encode(val).decode('unicode_escape')
    if typ.startswith('S'):    return str.encode(val).decode('unicode_escape').encode("ascii")
    return np.asscalar( np.array((val)).astype(typ) )

def check_destination_exists(input_file, output_dir):
	with pysam.AlignmentFile(input_file, "rc") as samfile:
		for read in samfile.fetch(until_eof=True):
			fast5_filename = read.get_tag(FILENAME_TAG)
			fast5_pathname = os.path.join(output_dir,fast5_filename)
			if os.path.exists(fast5_pathname):
				sys.exit( "Destination file already exists:{}".format(fast5_pathname) )

def is_empty_hdf_node(hdf_path):
	return hdf_path.endswith("/dummy_attr")
