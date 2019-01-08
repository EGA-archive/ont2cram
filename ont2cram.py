import os
import sys
import h5py
import pysam
import tqdm
import array
import argparse
import re
import numpy

FIRST_TAG    = "a0"
LAST_TAG     = "zZ"
FILENAME_TAG = "X0"
RESERVED_TAGS = [FILENAME_TAG]

class Tag:
    DIGITS = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    BASE = len(DIGITS)
    def __init__(self, start_tag="00"):
        self.current_tag_num = self.tag_to_int(start_tag)

    def tag_to_int(self, tag_name):
        return self.DIGITS.index(tag_name[0])*self.BASE+self.DIGITS.index(tag_name[1])

    def increment(self):
        self.current_tag_num += 1

    def get_name(self):
        return self.int_to_tag(self.current_tag_num)

    def int_to_tag(self, num):
        assert num >= 0
        base = len(self.DIGITS)
        res = ""
        while not res or num > 0:
            num, i = divmod(num, base)
            res = self.DIGITS[i] + res
        return res

global_dict_attributes = {}

def get_tag_type(np_typecode):
    if np_typecode.startswith(('u','i')): return 'i'
    if np_typecode.startswith( 'f'     ): return 'f'
    if np_typecode.startswith(('S','U')): return None
    sys.exit("Unknown type:{}".format(np_typecode))

def convert_t(typ):
    if typ.startswith(('=','>','<','|')): return typ[1:]
    return typ

def get_type(value):
    try:
        return value.dtype
    except AttributeError:
        return  numpy.dtype(type(value))

def convert_type(value):
    typ = convert_t( get_type(value).str ) 
    return ( value.decode('ascii') if typ[0]=='S' else value, typ )

def is_fastq_path(hdf_path):
    return hdf_path.endswith("BaseCalled_template/Fastq")
    
def is_signal_path(hdf_path):
    return "/Reads/Read" in hdf_path and hdf_path.endswith("Signal") 

def types_equal(t1,t2):
    if t1.startswith('S') and t2.startswith('S'): return True
    if t1.startswith('U') and t2.startswith('U'): return True
    return t1==t2 

def process_dataset(hdf_path, columns):
    if is_fastq_path(hdf_path): 
        return
    for column in columns:
        col_name = column[0]
        col_type_str = str(column[1][0]) if isinstance(column[1], tuple) else column[1]
        col_type = convert_t(col_type_str)
        full_key = hdf_path+'/'+col_name
        try:
            if not types_equal( global_dict_attributes[full_key][0], col_type ):
                sys.exit("Column '{}' - different types in fast5 files: {} vs {}".format(full_key,global_dict_attributes[full_key][0],col_type) )
        except KeyError:
            global_dict_attributes[full_key] = [col_type, 0]  
        

def remove_read_number(attribute_path):
    if "/Reads/Read_" in attribute_path: 
        return re.sub(r"(^.*Read_)(\d+)(.*$)", r"\g<1>XXX\g<3>", attribute_path)
    else:
        return attribute_path

    
def pre_process_group_attrs(node_path, hdf_node):
    global global_dict_attributes

    node_path = remove_read_number( node_path )    

    if type(hdf_node) is h5py.Dataset:
        columns = hdf_node.dtype.fields.items() if hdf_node.dtype.fields else [('noname', hdf_node.dtype.str)]
        process_dataset( node_path, columns )
           
    for key, val in hdf_node.attrs.items():
        full_key = node_path+'/'+key
        try:
            pair = global_dict_attributes[full_key]
            if pair[0] == val: pair[1] += 1
            global_dict_attributes[full_key] = pair
        except KeyError:
            global_dict_attributes[full_key] = [ val, 1 ]  

def walk_fast5( filename, walk_group_function ):
    with h5py.File(filename,'r') as f:
        walk_group_function("/", f)
        f.visititems( walk_group_function )

def get_list_of_fast5_files( dir ):
    return [os.path.join(dir, f) for f in os.listdir(dir) if os.path.isfile(os.path.join(dir,f)) and f.endswith('.fast5')]


def is_shared_value(value, total_fast5_files):
    return value > total_fast5_files//2


def write_cram(fast5_files, cram_file, skipsignal):
    total_fast5_files = len(fast5_files)
    comments_list = []
    tag = Tag(FIRST_TAG)
    for key,val in global_dict_attributes.items():

        is_column = True if val[1]==0 else False

        value,hdf_type = (None,val[0]) if is_column else convert_type(val[0])

        tag_and_val = "TG:"+tag.get_name()
        tag.increment()
        if tag_and_val.endswith(LAST_TAG): sys.exit("Running out of Tag space : too many atributes in Fast5")
        
        if is_shared_value(val[1], total_fast5_files) and "read_number" not in key:  
            tag_and_val += " CV:"+repr(value)

        comments_list.append( "{}:'{}':{} {}".format( "COL" if is_column else "ATR", key, hdf_type, tag_and_val) )
        global_dict_attributes[key][1] = tag_and_val

            
    header = {  'HD': {'VN': '1.0'},
                'SQ': [{'LN': 0, 'SN': '*'}],
                'CO': comments_list 
               }

    with pysam.AlignmentFile( cram_file, "wc", header=header, format_options=[b"no_ref=1"] ) as outf:
        for filename in tqdm.tqdm(fast5_files): 
            with h5py.File(filename,'r') as fast5:
                a = pysam.AlignedSegment()

                def get_tag_name_cv(hdf_full_path):
                    pair = global_dict_attributes[hdf_full_path]
                    assert( pair[1].startswith("TG:") )
                    tag_name = pair[1][3:5]
                    pos_CV = pair[1].find("CV:")
                    val_CV = None if pos_CV==-1 else pair[1][pos_CV+3:]  
                    return (tag_name,val_CV)              

                def get_column( dset, col_name ):
                    if col_name=="noname": return dset[()]
                    return dset[col_name]
                    
                def process_dataset(hdf_path, dset, columns):
                    if is_signal_path(hdf_path) and skipsignal: return
                    if is_fastq_path(hdf_path)                : return
                    for column in columns:
                        col_name   = column[0]
                        tag_name,_ = get_tag_name_cv(hdf_path+'/'+col_name)
                        col_array  = get_column(dset,col_name).tolist()
                        a.set_tag(tag_name, b''.join(col_array) if type(col_array[0]) is bytes else col_array)
                                        
                fastq_path  = None
                def process_attrs( name, group_or_dset ):
                    nonlocal fastq_path
                    if is_fastq_path(name): fastq_path=name
                    name = remove_read_number( name )    

                    if type(group_or_dset) is h5py.Dataset:
                        columns = group_or_dset.dtype.fields.items() if group_or_dset.dtype.fields else [('noname', None)]
                        process_dataset( name, group_or_dset, columns )
                    
                    for key, val in group_or_dset.attrs.items():
                        value, hdf_type = convert_type(val)
                        tag_name, val_CV = get_tag_name_cv(name+'/'+key)
                        try:
                            if repr(value) != val_CV : a.set_tag( tag_name, value, get_tag_type(hdf_type) )
                        except ValueError:
                            sys.exit("Could not detemine tag type (val={}, hdf_type={})".format(value,hdf_type))
                                                            
                fast5.visititems( process_attrs )

                if not fastq_path: 
                    sys.exit("Bad Fast5: Fastq dataset could not be found in '{}'".format(filename))

                a.set_tag( FILENAME_TAG, os.path.basename(filename) )
                      
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
    parser.add_argument('-s','--skipsignal', help='Skips the raw signal data', action='store_true')
    args = parser.parse_args()

    if not os.path.isdir(args.inputdir): sys.exit( 'Not dir: %s' % args.inputdir )

    fast5_files =  get_list_of_fast5_files( args.inputdir )

    print("Phase 1 of 2 : pre-processing Fast5 files...")
    for f in tqdm.tqdm(fast5_files): 
        walk_fast5( f, pre_process_group_attrs )

    print("Phase 2 of 2 : converting Fast5 files to CRAM..." )
    write_cram( fast5_files, args.outputfile, args.skipsignal )

if __name__ == "__main__":
    main()