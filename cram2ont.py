import os
import re
import sys
import h5py
import tqdm
import shlex
import pysam
import numpy
import argparse
import numpy.lib.recfunctions as rfn

DT_STR_VLEN = h5py.special_dtype(vlen=str)

FILENAME_TAG = "X0"
RESERVED_TAGS = [FILENAME_TAG]
STR_HEX_PATTERN = re.compile(r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$")

def convert_t(typ):
    if typ=='f16':  return numpy.float16
    if typ=='f32':  return numpy.float32
    if typ=='f64':  return numpy.float64
    if typ=='i8':   return numpy.int8
    if typ=='i16':  return numpy.int16
    if typ=='i32':	return numpy.int32
    if typ=='i64':	return numpy.int64
    if typ=='u8':	return numpy.uint8
    if typ=='u16':	return numpy.uint16
    if typ=='u32':	return numpy.uint32
    if typ=='u64':	return numpy.uint64
    if typ.startswith('U'):  return typ
    if typ.startswith('S'):  return typ

def convert_type(val, typ):
    if typ.startswith('U'):    return str.encode(val).decode('unicode_escape').encode("utf-8")
    if typ.startswith('S'):    return str.encode(val).decode('unicode_escape').encode("ascii")
    return convert_t(typ)(val)

def join_struct_arrays(arrays):
    if len(arrays)==1: return arrays[0]
    newdtype = sum((a.dtype.descr for a in arrays), [])
    newrecarray = numpy.empty(len(arrays[0]), dtype = newdtype)
    for a in arrays:
        for name in a.dtype.names:
            newrecarray[name] = a[name]
    return newrecarray
    
def cram_to_fast5(cram_filename, output_dir):
    class Attribute:
      path = ''
      type = ''
      value = None
      is_col = False

    attr_dict = {}    
    with pysam.AlignmentFile(cram_filename, "rc") as samfile:
        read_number_tag = None
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

            if a.path == "Raw/Reads/Read_XXX/read_number":
                read_number_tag = tag
        
        if not read_number_tag:
            sys.exit("Could not find read_number")

        def is_hex_str( obj, expected_len, is_null_term=False ):
            bytes_to_remove_end = 1 if is_null_term else 0
            return isinstance(obj,bytes) and len(obj)==expected_len and (not is_null_term or obj[-1:]==0) \
                   and STR_HEX_PATTERN.match(obj[:expected_len-bytes_to_remove_end].decode())

        def get_path(hdf_path, read_number):
            return hdf_path.replace("Read_XXX",  read_number)       
            
        def write_hdf_attr(hdf5_file, full_attr_path, attr_value):
            group_name,_,attr_name =  full_attr_path.rpartition('/') 
            if attr_name=="noname": raise
            try:
                group = hdf5_file[group_name]
            except KeyError:               
                group = hdf5_file.create_group(group_name)

            if is_hex_str(attr_value,36):
                group.attrs.create(attr_name, attr_value, dtype="|S36")
            else:
                group.attrs[attr_name] = attr_value 
            
        for read in tqdm.tqdm(samfile.fetch(until_eof=True)):
            fast5_filename = read.get_tag(FILENAME_TAG)
            read_number  = "Read_"+str(read.get_tag(read_number_tag))
             
            with h5py.File( os.path.join(output_dir,fast5_filename), "w" ) as f:
                fastq_lines = numpy.string_(
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
                        #if col_name=="noname": print( "dset={}, col={}, type={}, val={}".format(dset_name,col_name,a.type,tag_val[:11]) )
                        #if a.type.startswith(('S','U')): print(tag_val)
                        #aa = list(chunkstring(tag_val,int(a.type[1:]))) if a.type.startswith(('S','U')) else tag_val
                        #print("fixed val=", aa[:11])
                                                
                        dset.append(
                            numpy.array( 
                                list(chunkstring(tag_val,int(a.type[1:]))) if a.type.startswith(('S','U')) else tag_val, 
                                dtype=None if col_name=="noname" else [(col_name, convert_t(a.type))] 
                                )
                             )
                for dset_name,columns in DSETS.items():
                    #f.create_dataset( dset_name, data=rfn.merge_arrays(columns, flatten=True, usemask=False) )
                    #print(f"Saving dset={dset_name}, cols={columns}, len={len(columns)} \n\n") 
                    a = rfn.merge_arrays(columns, flatten=True, usemask=False)
                    f.create_dataset( dset_name, data=a )
                                        
                # write constant values stored in cram header
                for a in attr_dict.values():
                    if a.is_col: continue
                    if a.value: write_hdf_attr( f, get_path(a.path,read_number), convert_type(a.value,a.type) )                     

                # write tags stored in cram records                                                
                for tag_name, tag_val in read.get_tags():
                    if tag_name in RESERVED_TAGS: continue    
                    a = attr_dict[tag_name]
                    if a.is_col: continue
                    #print( "or_t={}, type={}, val={}".format(a.type, str(type(tag_val)), tag_val) )
                    if a.value != tag_val: 
                        write_hdf_attr( f, get_path(a.path,read_number), convert_type(tag_val,a.type) )


def main():
    parser = argparse.ArgumentParser(description='CRAM to Fast5 conversion utility (this is a reverse converter which allows to restore original Fast5 collection from Cram generated by ont2cram)')
    parser.add_argument( '-i','--inputfile', help='Input CRAM filename', required=True )
    parser.add_argument( '-o','--outputdir', help='Output directory for generated Fast5 files', action='store', default=os.getcwd() )
    args = parser.parse_args()

    if not os.path.isdir(args.outputdir): 
        sys.exit( 'Not dir: %s' % args.outputdir )

    cram_to_fast5( args.inputfile, args.outputdir )

if __name__ == "__main__":
    main()