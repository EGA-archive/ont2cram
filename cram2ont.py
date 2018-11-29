import os
import sys
import h5py
import shlex
import pysam
import argparse

FILENAME_TAG = "X0"
SIGNAL_TAG   = "zZ"
RESERVED_TAGS = [SIGNAL_TAG, FILENAME_TAG]

def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))
    
def cram_to_fast5(cram_filename, output_dir):
    class Attribute:
      path = ''
      type = ''
      value = None
      
    attr_dict = {}
    with pysam.AlignmentFile(cram_filename, "rc") as samfile:
        filename_tag = None
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

        def write_hdf_attr(hdf5_file, full_attr_path, attr_value):
            group_name,_,attr_name =  full_attr_path.rpartition('/')                
            group = hdf5_file.require_group(group_name)
            group.attrs[attr_name] = attr_value
            
        def encode_attr(value):
            return value.encode("ascii") if a.type=="bytes" else value
            
        for read in samfile.fetch(until_eof=True):
            fast5_filename = read.get_tag(FILENAME_TAG)
             
            with h5py.File( os.path.join(output_dir,fast5_filename), "w" ) as f:
                for a in attr_dict.values():
                    if a.value: write_hdf_attr( f, a.path, encode_attr(a.value) ) 

                for tag_name, tag_val in read.get_tags():
                    if tag_name in RESERVED_TAGS: continue    
                    a = attr_dict[tag_name]
                    #print( "or_t={}, type={}, val={}".format(a.type, str(type(tag_val)), tag_val) )

                    if a.value != tag_val: 
                        write_hdf_attr( f, a.path, encode_attr(tag_val) )
            #print(read)


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