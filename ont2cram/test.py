#!/usr/bin/env python3
import os
import sys
import shutil
import tempfile
import unittest
import subprocess
import argparse
from parameterized import parameterized

from ont2cram.converter import converter
from ont2cram.reverse_converter import reverse_converter

KEEP_TMP = False
IGNORE_LINES  = [
	'.fast5" {',
	"H5T_VARIABLE",
	"STRPAD",
	"H5T_CSET_UTF8",
	"DATASPACE"
]

TEST_DATA_DIR = os.path.join(os.getcwd(),"test_data") ###################

def test(keep_tmp=False):

	global KEEP_TMP
	KEEP_TMP = keep_tmp

	del(sys.argv[1:])
	unittest.main(verbosity=2, buffer=not KEEP_TMP)

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
    	self.fast5_original_dir   = TEST_DATA_DIR

    def setUp(self):
        self.fast5_restored_dir  = tempfile.mkdtemp()
        self.h5dump_original_dir = tempfile.mkdtemp()
        self.h5dump_restored_dir = tempfile.mkdtemp()

    def tearDown(self):
    	if KEEP_TMP:
    		print("----------------------TMP DIRS-----------------------")
    		print("Fast5  original dir:'{}'".format(self.fast5_original_dir))
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
                error_msg = "File:'{}'".format(os.path.basename(path1))

                if KEEP_TMP:
                	error_msg = 'diff -u "{}" "{}" | cdiff'.format(path1,path2)
                	print(error_msg)

                for line1,line2,n in zip(f1,f2,range(sys.maxsize)):
                    if any(x in line2 for x in IGNORE_LINES):
                    	continue
                    buf1.append(line1)
                    buf2.append(line2)
                    if n%10000==0:
                        self.assertEqual(buf1,buf2, error_msg )
                        buf1.clear()
                        buf2.clear()

                self.assertEqual(buf1,buf2, error_msg )
                buf1.clear()
                buf2.clear()


    def assert2DirsEqual(self, dir1, dir2):
        dir1_files = os.listdir(dir1)
        dir2_files = os.listdir(dir2)
        self.assertTrue( len(dir1)==len(dir2) )
        if KEEP_TMP: print("---------------------DIFF FILES----------------------")
        for f in dir1_files:
            self.assert2FilesEqual( os.path.join(dir1,f), os.path.join(dir2,f) )

    @parameterized.expand([x[0] for x in os.walk(TEST_DATA_DIR)])
    def test_round_trip(self, fast5_dir):
        self.fast5_original_dir = fast5_dir
        files = get_all_fast5_files(self.fast5_original_dir)
        print("Testing dir={}, total files={}".format(self.fast5_original_dir, len(files)))
        if len(files)==0: return

        cram_path = tempfile.mkstemp(suffix=".cram")[1]
        try:

            #forward conversion
            converter(self.fast5_original_dir, None, cram_path, False )

            #reverse conversion
            reverse_converter( cram_path, self.fast5_restored_dir )

            h5dump_all_files_in_dir( self.fast5_original_dir,  self.h5dump_original_dir )
            h5dump_all_files_in_dir( self.fast5_restored_dir, self.h5dump_restored_dir )

            self.assert2DirsEqual( self.h5dump_original_dir, self.h5dump_restored_dir )

        finally:
        	if KEEP_TMP:
        		print("----------------------TMP FILES----------------------")
        		print("Generated CRAM :'{}'".format(cram_path))
        	else:
        		os.remove(cram_path)
