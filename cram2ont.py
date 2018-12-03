import os
import re
import sys
import h5py
import tqdm
import shlex
import pysam
import numpy
import argparse

FILENAME_TAG = "X0"
SIGNAL_TAG   = "zZ"
RESERVED_TAGS = [SIGNAL_TAG, FILENAME_TAG]
STR_HEX_PATTERN = re.compile(r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$")

def convert_type(val, typ):
    if typ=='f32':  return numpy.float32(val)
    if typ=='f64':  return numpy.float64(val)
    if typ=='i8':   return numpy.int8(val)
    if typ=='i32':	return numpy.int32(val)
    if typ=='i64':	return numpy.int64(val)
    if typ=='ui8':	return numpy.uint8(val)
    if typ=='ui32':	return numpy.uint32(val)
    if typ=='ui64':	return numpy.uint64(val)
    if typ=='str':  return str.encode(val).decode('unicode_escape').encode("utf-8")
    if typ=='bytes':return str.encode(val).decode('unicode_escape').encode("ascii")
    
def cram_to_fast5(cram_filename, output_dir):
    class Attribute:
      path = ''
      type = ''
      value = None

    attr_dict = {}    
    with pysam.AlignmentFile(cram_filename, "rc") as samfile:
        read_number_tag = None
        for comment in samfile.header["CO"]:
            if not comment.startswith("ONT:"): continue

            parts = shlex.split(comment)
            part0 = parts[0].split(':')

            a=Attribute()
            a.path = part0[1]
            a.type = part0[2] 
            a.value = parts[2][3:] if len(parts)==3 else None

            tag = parts[1][3:]
            
            attr_dict[tag]=a

            if a.path == "Raw/Reads/Read_XXX/read_number":
                read_number_tag = tag
        
        if not read_number_tag:
            sys.exit("Could not find read_number")

        def is_hex_str( obj, expected_len, is_null_term=False ):
            bytes_to_remove_end = 1 if is_null_term else 0
            return isinstance(obj,bytes) and len(obj)==expected_len and (not is_null_term or obj[-1:]==0) \
                   and STR_HEX_PATTERN.match(obj[:expected_len-bytes_to_remove_end].decode())           
            
        def write_hdf_attr(hdf5_file, full_attr_path, attr_value, read_number):
            full_attr_path = full_attr_path.replace("Read_XXX",  read_number)
            group_name,_,attr_name =  full_attr_path.rpartition('/')                
            group = hdf5_file.require_group(group_name)

            if is_hex_str(attr_value,36):
                group.attrs.create(attr_name, attr_value, dtype="|S36")
            else:
                group.attrs[attr_name] = attr_value 
            
        for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
            fast5_filename = read.get_tag(FILENAME_TAG)
            read_number  = "Read_"+str(read.get_tag(read_number_tag))
             
            with h5py.File( os.path.join(output_dir,fast5_filename), "w" ) as f:
                for a in attr_dict.values():
                    if a.value: write_hdf_attr( f, a.path, convert_type(a.value,a.type), read_number ) 

                for tag_name, tag_val in read.get_tags():
                    if tag_name in RESERVED_TAGS: continue    
                    a = attr_dict[tag_name]
                    #print( "or_t={}, type={}, val={}".format(a.type, str(type(tag_val)), tag_val) )

                    if a.value != tag_val: 
                        write_hdf_attr( f, a.path, convert_type(tag_val,a.type), read_number )

def main():
    parser = argparse.ArgumentParser(description='CRAM to Fast5 conversion utility')
    parser.add_argument( '-i','--inputfile', help='Input CRAM filename', required=True )
    parser.add_argument( '-o','--outputdir', help='Output directory for generated Fast5 files', action='store', default=os.getcwd() )
    args = parser.parse_args()

    if not os.path.isdir(args.outputdir): 
        sys.exit( 'Not dir: %s' % args.outputdir )

    cram_to_fast5( args.inputfile, args.outputdir )

if __name__ == "__main__":
    main()