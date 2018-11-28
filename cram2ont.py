import os
import sys
import h5py
import shlex
import pysam
import argparse

def dump(obj):
  for attr in dir(obj):
    print("obj.%s = %r" % (attr, getattr(obj, attr)))
    
def cram_to_fast5( cram_filename, output_dir):
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

            if a.path == "UniqueGlobalKey/context_tags/filename":
                filename_tag = tag 

        # for x in attr_dict.values():
        #    print(x.value) 

        if not filename_tag: 
            sys.exit("Could not detect original Fast5 filename")

        def write_hdf_attr(hdf5_file, full_attr_path, attr_value):
            group_name,_,attr_name =  full_attr_path.rpartition('/')                
            group = hdf5_file.require_group(group_name)
            group.attrs[attr_name] = attr_value
            
        for read in samfile.fetch(until_eof=True):
            fast5_filename = read.get_tag(filename_tag)
             
            with h5py.File(fast5_filename, "w") as f:
                for a in attr_dict.values():
                    if a.value: write_hdf_attr(f, a.path, a.value)

                for tag_name, tag_val in read.get_tags():
                    if tag_name == "zz": continue    
                    a = attr_dict[tag_name]
                    if a.value != tag_val: write_hdf_attr(f, a.path, tag_val)
            #print(read)


def main():
    parser = argparse.ArgumentParser(description='CRAM to Fast5 conversion utility')
    parser.add_argument( '-i','--inputfile', help='Input CRAM filename', required=True )
    parser.add_argument( '-o','--outputdir', help='Output directory for generated Fast5 files', action='store', default=os.getcwd() )
    args = parser.parse_args()

    if not os.path.isdir(args.outputdir): sys.exit( 'Not dir: %s' % args.inputdir )

    cram_to_fast5( args.inputfile, args.outputdir )

    #print("Phase 2 of 2 : converting Fast5 files to CRAM..." )

if __name__ == "__main__":
    main()