import ont2cram
import cram2ont

import unittest
import os, subprocess
import shutil, tempfile

def h5dump_all_files_in_dir(input_dir, output_dir):
    for f in os.listdir(input_dir):
        file_path = os.path.join(input_dir ,f)
        dump_path = os.path.join(output_dir,f+".txt") 
    with open(dump_path,'w') as outfile:
        subprocess.call(["h5dump",file_path], stdout=outfile)

def assert2FilesEqual(path1, path2):
    with open(path1) as f1:
        with open(path2) as f2:    
            assert( [row for row in f1] == [row for row in f2] )
    
def assert2DirsEqual(dir1, dir2):
    dir1_files = os.listdir(dir1)
    dir2_files = os.listdir(dir2)
    assert( len(dir1)==len(dir2) )
    for f in dir1_files: assert2FilesEqual( os.path.join(dir1,f), os.path.join(dir2,f) )      

class Ont2CramTests(unittest.TestCase):
    def setUp(self):
        self.fast5_origial_dir   = os.path.join(os.getcwd(),"test")
        self.fast5_restored_dir  = tempfile.mkdtemp()        
        self.h5dump_original_dir = tempfile.mkdtemp()        
        self.h5dump_restored_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fast5_restored_dir)
        shutil.rmtree(self.h5dump_original_dir)
        shutil.rmtree(self.h5dump_restored_dir)                                                                

    def test_round_trip(self):
        cram_path = tempfile.mkstemp()[1]        
        try:
            #forward conversion
            ont2cram.main( self.fast5_origial_dir, None, cram_path, False ) 
            #reverse conversion
            cram2ont.main( cram_path, self.fast5_restored_dir )
            
            h5dump_all_files_in_dir( self.fast5_origial_dir,  self.h5dump_original_dir )    
            h5dump_all_files_in_dir( self.fast5_restored_dir, self.h5dump_restored_dir )   

            assert2DirsEqual( self.h5dump_original_dir, self.h5dump_restored_dir )            
        finally:
            os.remove(cram_path)           
