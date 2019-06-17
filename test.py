#!/usr/bin/env python3
import ont2cram
import cram2ont

import os
import sys
import shutil
import tempfile
import unittest
import subprocess
import argparse

KEEP_TMP      = False
IGNORE_LINES  = [ 
	"HDF5", 
	"STRSIZE", 
	"STRPAD", 
#	"ASCII", 
#	"H5T_STD_I16LE", 
	"DATASPACE"
]

TEST_DATA_DIR = os.path.join(os.getcwd(),"test_data")

def get_all_fast5_files(dir):
	return  [f for f in os.listdir(dir) if os.path.isfile(os.path.join(dir, f)) and f.endswith(".fast5")]	

def h5dump_all_files_in_dir(input_dir, output_dir):
    files = get_all_fast5_files(input_dir)
    for f in files:
        file_path = os.path.join(input_dir ,f)
        dump_path = os.path.join(output_dir,f+".txt")
        with open(dump_path,'w') as outfile:
    		# Cram does not support double precision floats, 
    		# so we limit here the dumped precision using the format string
        	subprocess.call(["h5dump","-m%.4g",file_path], stdout=outfile)
    

class Ont2CramTests(unittest.TestCase):

    def __init__(self, *args, **kwargs):
    	super( Ont2CramTests, self ).__init__(*args, **kwargs)
    	self.fast5_origial_dir   = TEST_DATA_DIR
    		
    def setUp(self):        
        self.fast5_restored_dir  = tempfile.mkdtemp()        
        self.h5dump_original_dir = tempfile.mkdtemp()        
        self.h5dump_restored_dir = tempfile.mkdtemp()

    def tearDown(self):
    	if KEEP_TMP:
    		print("-----TMP DIRS-----")
    		print("Fast5  original dir:'{}'".format(self.fast5_origial_dir))
    		print("Fast5  restored dir:'{}'".format(self.fast5_restored_dir))
    		print("h5dump original dir:'{}'".format(self.h5dump_original_dir))
    		print("h5dump restored dir:'{}'".format(self.h5dump_restored_dir))
    	else:
    		shutil.rmtree(self.fast5_restored_dir)
    		shutil.rmtree(self.h5dump_original_dir)
    		shutil.rmtree(self.h5dump_restored_dir)                                                                

    def assert2FilesEqual(self, path1, path2):
        with open(path1) as f1:
            with open(path2) as f2:
                buf1 = []
                buf2 = []        
                n=0    
                error_msg = "'{}' vs'{}' (line:{})"
                for line1,line2  in zip(f1,f2):
                    n += 1
                    if any(x in line1 for x in IGNORE_LINES): 
                    	continue
                    buf1.append(line1.replace("Read_","read_"))
                    buf2.append(line2)                    
                    if n%10000==0:                    
                        self.assertEqual(buf1,buf2, error_msg.format(path1,path2,n))                        
                        buf1.clear()
                        buf2.clear()
                self.assertEqual(buf1,buf2, error_msg.format(path1,path2,n))                        
                                
                   
    def assert2DirsEqual(self, dir1, dir2):
        dir1_files = os.listdir(dir1)
        dir2_files = os.listdir(dir2)
        self.assertTrue( len(dir1)==len(dir2) )
        for f in dir1_files: 
            self.assert2FilesEqual( os.path.join(dir1,f), os.path.join(dir2,f) )                 

    def round_trip_test(self, fast5_dir):
        self.fast5_origial_dir = fast5_dir
        files = get_all_fast5_files(self.fast5_origial_dir)
        print("Testing dir={}, files={}".format(self.fast5_origial_dir,files))        
        if len(files)==0: return
                
        cram_path = tempfile.mkstemp()[1]        
        try:

            #forward conversion
            ont2cram.run( self.fast5_origial_dir, None, cram_path, False ) 
            
            #reverse conversion
            cram2ont.run( cram_path, self.fast5_restored_dir )
            
            h5dump_all_files_in_dir( self.fast5_origial_dir,  self.h5dump_original_dir )    
            h5dump_all_files_in_dir( self.fast5_restored_dir, self.h5dump_restored_dir )   

            self.assert2DirsEqual( self.h5dump_original_dir, self.h5dump_restored_dir )            
        finally:
            os.remove(cram_path)

    def test_all(self):    	
    	subdirs = [x[0] for x in os.walk(self.fast5_origial_dir)]
    	for dir in subdirs:
    		with self.subTest():
    			self.round_trip_test(dir)    		

def main():
	parser = argparse.ArgumentParser(description='Fast5 to CRAM testing tool')
	parser.add_argument('-k','--keeptmp', help='Do not delete tmp folders and files on exit', action='store_true')
	args = parser.parse_args()

	global KEEP_TMP
	KEEP_TMP = args.keeptmp	
	
	del(sys.argv[1:])
	unittest.main()

if __name__ == '__main__':
	main()