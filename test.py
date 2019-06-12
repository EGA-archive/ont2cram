#!/usr/bin/env python3
import ont2cram
import cram2ont

import os
import sys
import shutil
import tempfile
import unittest
import subprocess

IGNORE_LINES  = ["HDF5"]#,"\n"]
TEST_DATA_DIR = os.path.join(os.getcwd(),"test_data")

def h5dump_all_files_in_dir(input_dir, output_dir):
    for f in os.listdir(input_dir):
        file_path = os.path.join(input_dir ,f)
        dump_path = os.path.join(output_dir,f+".txt") 
    with open(dump_path,'w') as outfile:
        subprocess.call(["h5dump",file_path], stdout=outfile)
    

class Ont2CramTests(unittest.TestCase):
    def setUp(self):
        self.fast5_origial_dir   = TEST_DATA_DIR
        self.fast5_restored_dir  = tempfile.mkdtemp()        
        self.h5dump_original_dir = tempfile.mkdtemp()        
        self.h5dump_restored_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fast5_restored_dir)
        shutil.rmtree(self.h5dump_original_dir)
        shutil.rmtree(self.h5dump_restored_dir)                                                                

    def assert2FilesEqual(self, path1, path2):
        #self.assertTrue(False,"rrrrrrr")
        
        with open(path1) as f1:
            with open(path2) as f2:
                buf1 = []
                buf2 = []        
                n=0    
                for line1,line2  in zip(f1,f2):
                    n += 1
                    if any(x in line1 for x in IGNORE_LINES): continue
                    if n%10==0:                    
                        #diff = d.compare("".join(buf1), "".join(buf2))
                        #print ''.join(list(diff))
                        self.assertEqual(buf1,buf2, path1)                        
                        buf1.clear()
                        buf2.clear()
                    else:
                        buf1.append(line1)
                        buf2.append(line2)
                self.assertEqual(buf1,buf2)                        
                   
    def assert2DirsEqual(self, dir1, dir2):
        dir1_files = os.listdir(dir1)
        dir2_files = os.listdir(dir2)
        self.assertTrue( len(dir1)==len(dir2) )
        for f in dir1_files: 
            self.assert2FilesEqual( os.path.join(dir1,f), os.path.join(dir2,f) )                 

    def test_round_trip(self):
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

def main():
    #del(sys.argv[1:])
    unittest.main()

if __name__ == '__main__':
    main()