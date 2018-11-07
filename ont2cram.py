import os
import sys
import h5py
import pysam
import tqdm
import array
import argparse
import numpy

global_dict_attributes = {}

def convert_type(value):
    typ= {
        "<class 'numpy.float32'>" : 'f32',
        "<class 'numpy.float64'>" : 'f64',                
        "<class 'numpy.int8'>"    : 'i8',
        "<class 'numpy.int32'>"   : 'i32',
        "<class 'numpy.int64'>"   : 'i64',  
        "<class 'numpy.uint8'>"   : 'ui8',
        "<class 'numpy.uint32'>"  : 'ui32',        
        "<class 'numpy.uint64'>"  : 'ui64',
        "<class 'str'>"           : 'str',
        "<class 'numpy.bytes_'>"  : 'bytes',
        "<class 'bytes'>"         : 'bytes'
    }[str(type(value))]
    return ( value.decode('ascii') if typ=='bytes' else value, typ )
    
def pre_process_group_attrs(name, group):
    global global_dict_attributes
    if "/Reads/Read" in name: return    
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
    tag =  int('a2', 36)
    for key,val in global_dict_attributes.items():
        (value, typ) = convert_type(val[0])
        tag_or_val = "CV:"+repr(value) 
        if not is_shared_value(val[1], total_fast5_files): 
            tag_or_val = "TG:"+numpy.base_repr(tag, 36).lower()
            tag += 1 

        comments_list.append( "ONT:'{}':{} {}".format(key, typ, tag_or_val) )
        global_dict_attributes[key][1] = tag_or_val
            
    header = {  'HD': {'VN': '1.0'},
                'SQ': [{'LN': 0, 'SN': '*'}],
                'CO': comments_list 
               }

    with pysam.AlignmentFile( cram_file, "wc", header=header, format_options=[b"no_ref=1"] ) as outf:
        for filename in tqdm.tqdm(fast5_files): 
            with h5py.File(filename,'r') as fast5:   
                fastq_lines = fast5['Analyses/Basecall_1D_000/BaseCalled_template/Fastq'].value.splitlines()
                a = pysam.AlignedSegment()
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
                
                tagA0 = []
                signal_path = ""
                def process_attrs( name, group ):
                    nonlocal signal_path
                    if "/Reads/Read" in name and name.endswith("Signal"): 
                        signal_path=name
                    for key, val in group.attrs.items():
                        (value, typ) = convert_type(val)
                        full_key = name+'/'+key
                        try:
                            pair = global_dict_attributes[full_key]
                            if pair[1].startswith("TG:"): a.set_tag(pair[1][3:], value, 'i' if typ=="i64" else None)
                        except KeyError:
                            # have not seen this attr before - save at read level
                            tagA0.append( "{},{},{}".format(full_key, typ, value) )
                    
                                            
                fast5.visititems( process_attrs )
                a.set_tag( "a0",';'.join(tagA0) )
                a.set_tag( "a1", array.array('h',fast5[signal_path].value) )

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