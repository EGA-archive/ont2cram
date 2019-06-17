#!/usr/bin/env python3
import os
import re
import sys
import h5py
import tqdm
import shlex
import pysam
import argparse
import numpy as np
import numpy.lib.recfunctions as rfn

DT_STR_VLEN = h5py.special_dtype(vlen=str)

READ_NUM_TAG = "X0"
FILENAME_TAG = "X1"
RESERVED_TAGS = [READ_NUM_TAG, FILENAME_TAG]

STR_HEX_PATTERN = re.compile(r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$")

def convert_type(val, typ):
    #print(f"val={val}, typ={typ}, type={type(val)}")
    if typ.startswith('U'):    return str.encode(val).decode('unicode_escape')
    if typ.startswith('S'):    return str.encode(val).decode('unicode_escape').encode("ascii")
    return np.asscalar( np.array((val)).astype(typ) )

def check_destination_exists(cram_filename, output_dir):
	with pysam.AlignmentFile(cram_filename, "rc") as samfile:
		for read in samfile.fetch(until_eof=True):
			fast5_filename = read.get_tag(FILENAME_TAG)
			fast5_pathname = os.path.join(output_dir,fast5_filename)
			if os.path.exists(fast5_pathname):
				sys.exit( "Destination file already exists:{}".format(fast5_pathname) )

def cram_to_fast5(cram_filename, output_dir):
    check_destination_exists(cram_filename, output_dir)
    	
    class Attribute:
      path = ''
      type = ''
      value = None
      is_col = False

    attr_dict = {}    
    with pysam.AlignmentFile(cram_filename, "rc") as samfile:
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

        def get_path(hdf_path, read_number):
            return hdf_path.replace("Read_XXX",read_number).replace("read_XXX",read_number)       
            
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
            
        for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
            fast5_filename = read.get_tag(FILENAME_TAG)
            read_number  = "read_"+str(read.get_tag(READ_NUM_TAG))                        
            
            with h5py.File( os.path.join(output_dir,fast5_filename), "a" ) as f:
                if read.query_name != "nofastq":
                    fastq_lines = np.string_(
                        "\n".join( [read.query_name, read.query_sequence, '+', pysam.array_to_qualitystring(read.query_qualities)+'\n'] ) )
                    f.create_dataset( "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq", data=fastq_lines )

                def chunkstring(string, length):
                    return (string[0+i:length+i] for i in range(0, len(string), length))

                DSETS = {}
                for tag_name,tag_val in read.get_tags():
                    if tag_name in RESERVED_TAGS: continue 
                    a = attr_dict[tag_name]
                    if a.is_col:
                        dset_name,_,col_name =  get_path(a.path,read_number).rpartition('/')
                        if dset_name.endswith("Fastq") and col_name=="noname": continue
                        if dset_name not in DSETS: DSETS[dset_name] = []
                        dset = DSETS[dset_name]

                        if col_name=="noname":
                            #print(f"path={a.path}, val={tag_val[:11]}")
                            dset.append(tag_val)
                        else:       
                            dset.append(
                                np.array( 
                                    list(chunkstring(tag_val,int(a.type[1:]))) if a.type.startswith(('S','U')) else tag_val, 
                                    dtype=[(col_name, a.type)] 
                                )
                             )
                for dset_name,columns in DSETS.items():
                    d = columns[0] if len(columns)==1 else rfn.merge_arrays(columns, flatten=True, usemask=False)
                    f.create_dataset( dset_name, data=d )
                                        
                # write constant values stored in cram header
                for a in attr_dict.values():
                    if a.is_col: continue
                    if a.value : write_hdf_attr( f, get_path(a.path,read_number), a.value, a.type )                     

                # write tags stored in cram records                                                
                for tag_name, tag_val in read.get_tags():
                    if tag_name in RESERVED_TAGS: continue    
                    a = attr_dict[tag_name]
                    if a.is_col: continue
                    if a.value != tag_val: 
                        write_hdf_attr( f, get_path(a.path,read_number), tag_val, a.type )


def run(input_file, output_dir):
    if not os.path.isdir(output_dir): sys.exit( 'Not dir: %s' % output_dir )
    cram_to_fast5( input_file, output_dir )

def main():
	parser = argparse.ArgumentParser(description='CRAM to Fast5 conversion utility (this is a reverse converter which allows to restore original Fast5 collection from Cram generated by ont2cram)')
	parser.add_argument( '-i','--inputfile', help='Input CRAM filename', required=True )
	parser.add_argument( '-o','--outputdir', help='Output directory for generated Fast5 files', action='store', default=os.getcwd() )
	args = parser.parse_args()
	run( args.inputfile, args.outputdir )

if __name__ == "__main__":
    main()