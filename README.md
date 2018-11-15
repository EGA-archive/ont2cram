# ont2cram
Oxford Nanopore HDF/Fast5 to CRAM conversion tool

USAGE: 
python ont2cram.py -i INPUTDIR -o OUTPUTFILE

  -i INPUTDIR, --inputdir INPUTDIR
        
        Input directory containing Fast5 files
                        
  -o OUTPUTFILE, --outputfile OUTPUTFILE
                        
                        Output CRAM filename


Implementation details:

There is a mapping table in the header that maps ONT attributes to lowercase SAM aux tags eg:
ONT:'Analyses/Basecall_1D_000/Configuration/calibration_strand/genome_name':bytes TG:b5 CV:'Lambda_3.6kb'

general format is : ONT:'<hdf_attribute_pathname>':<original_datatype> TG:<2_letter_lowecase_tag> CV:<constant_value>

CV part is optional - currently present only if >50% of fast5 files have this value. Thus, in current implementation CV is more like 'common value' than constant ( it can be overwritten on the read level for those reads that have different value ).

Tag names are generated sequentially ( a0, a1, ...., aa, ab, ..... zz ). 
If 'zz' is reached the program exists with an error ( we need to discuss this problem since 2 letter tagspace is not always enougth to cover all HDF attributes )

Signal is stored in 'zz' tag

We need to discuss if we want to store Events(currently not stored) - James has mentioned that events are only present in old ONT data.

Some string values in HDF attributes have line breaks inside - is it valid for Cram tags or better to remove them?
