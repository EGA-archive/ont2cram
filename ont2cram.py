import os
import sys
import h5py
import pysam
import tqdm
import array
import argparse
import numpy
import re

global_dict_attributes = {}

def convert_type(value):
    typ= {
        "<class 'numpy.float32'>" : ('f32','f'),
        "<class 'numpy.float64'>" : ('f64','f'),               
        "<class 'numpy.int8'>"    : ('i8','i'),
        "<class 'numpy.int32'>"   : ('i32','i'),
        "<class 'numpy.int64'>"   : ('i64','i'),  
        "<class 'numpy.uint8'>"   : ('ui8','i'),
        "<class 'numpy.uint32'>"  : ('ui32','i'),        
        "<class 'numpy.uint64'>"  : ('ui64','i'),
        "<class 'str'>"           : ('str',None),
        "<class 'numpy.bytes_'>"  : ('bytes',None),
        "<class 'bytes'>"         : ('bytes',None)
    }[str(type(value))]
    return ( value.decode('ascii') if typ[0]=='bytes' else value, typ[0], typ[1] )

def remove_read_number(attribute_path):
    if "/Reads/Read_" in attribute_path: 
        return re.sub(r"(^.*Read_)(\d+)(.*$)", r"\g<1>XXX\g<3>", attribute_path)
    else:
        return attribute_path
    
def pre_process_group_attrs(name, group):
    global global_dict_attributes

    name = remove_read_number( name )    
           
    for key, val in group.attrs.items():
        full_key = name+'/'+key
        try:
            pair = global_dict_attributes[full_key]
            if( pair[0] == val ): pair[1] += 1
            global_dict_attributes[full_key] = pair
        except KeyError:
            global_dict_attributes[full_key] = [ val, 1 ]  

def walk_fast5( filename, walk_group_function ):
    with h5py.File(filename,'r') as f:
        f.visititems( walk_group_function )

def get_list_of_fast5_files( dir ):
    return [os.path.join(dir, f) for f in os.listdir(dir) if os.path.isfile(os.path.join(dir,f)) and f.endswith('.fast5')]


def is_shared_value(value, total_fast5_files):
    return value > total_fast5_files//2


def write_cram(fast5_files, cram_file):
    total_fast5_files = len(fast5_files)
    comments_list = []
    tag = int('a0', 36)
    for key,val in global_dict_attributes.items():
        value,hdf_type,_ = convert_type(val[0])

        tag_and_val = "TG:"+numpy.base_repr(tag, 36).lower()
        tag += 1 
        
        if is_shared_value(val[1], total_fast5_files):  tag_and_val += " CV:"+repr(value)

        comments_list.append( "ONT:'{}':{} {}".format(key, hdf_type, tag_and_val) )
        global_dict_attributes[key][1] = tag_and_val

        if tag == "zz" : sys.exit("Running out of Tag space : too many atributes in Fast5")
            
    header = {  'HD': {'VN': '1.0'},
                'SQ': [{'LN': 0, 'SN': '*'}],
                'CO': comments_list 
               }

    with pysam.AlignmentFile( cram_file, "wc", header=header, format_options=[b"no_ref=1"] ) as outf:
        for filename in tqdm.tqdm(fast5_files): 
            with h5py.File(filename,'r') as fast5:
                a = pysam.AlignedSegment()

                signal_path = None
                fastq_path  = None

                def process_attrs( name, group ):
                
                    nonlocal signal_path
                    nonlocal fastq_path
                    if "/Reads/Read" in name and name.endswith("Signal"): signal_path=name
                    if name.endswith("BaseCalled_template/Fastq")       : fastq_path=name

                    name = remove_read_number( name )    
                    
                    for key, val in group.attrs.items():
                        (value, hdf_type, tag_type) = convert_type(val)
                        full_key = name+'/'+key
                        pair = global_dict_attributes[full_key]
                        assert( pair[1].startswith("TG:") )

                        tag_name = pair[1][3:5]
                        pos_CV = pair[1].find("CV:")
                        val_CV = None if pos_CV==-1 else pair[1][pos_CV+3:]
                        try:
                            if repr(value) != val_CV : a.set_tag(tag_name, value, tag_type)
                        except ValueError:
                            sys.exit("Could not detemine tag type (val={}, hdf_type={})".format(value,hdf_type))
                                                            
                fast5.visititems( process_attrs )

                if( not signal_path or not fastq_path ): 
                    sys.exit("Bad Fast5: signal or fastq could not be found in '{}'".format(filename))

                a.set_tag( "zz", array.array('h',fast5[signal_path].value) )
                      
                fastq_lines = fast5[fastq_path].value.splitlines()
                
                a.query_name = fastq_lines[0]
                a.query_sequence=fastq_lines[1]
                a.flag = 4
                a.reference_id = 0
                a.reference_start = 0
                a.mapping_quality = 0
                a.cigar = ()
                a.next_reference_id = 0
                a.next_reference_start=0
                a.template_length=0
                a.query_qualities = pysam.qualitystring_to_array(fastq_lines[3])

                outf.write(a)

def main():
    parser = argparse.ArgumentParser(description='Fast5 to CRAM conversion utility')
    parser.add_argument('-i','--inputdir', help='Input directory containing Fast5 files', required=True)
    parser.add_argument('-o','--outputfile', help='Output CRAM filename', required=True)
    args = parser.parse_args()

    if not os.path.isdir(args.inputdir): sys.exit( 'Not dir: %s' % args.inputdir )

    fast5_files =  get_list_of_fast5_files( args.inputdir )

    print("Phase 1 of 2 : pre-processing Fast5 files...")
    for f in tqdm.tqdm(fast5_files): 
        walk_fast5( f, pre_process_group_attrs )

    print("Phase 2 of 2 : converting Fast5 files to CRAM..." )
    write_cram( fast5_files, args.outputfile )

if __name__ == "__main__":
    main()