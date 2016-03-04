uPEP-helpers
============

A core fabfile for managing uPEP databases

Mitchell Stanton-Cook (m.stantoncook@gmail.com)

Requirements (python libraries):
    * ftputil
    * fabric
    * mysql.connector

Requirements (3rd party):
    * blast2 (legacy)
    * MySQL

Overview
--------

Aids in updating uPEP web app databases. If a new RefSeq release will:
    * download all required NCBI RefSeq databases, 
    * compact the databases, 
    * compile the databases, 
    * apply uPEP finding (creating uPEP databases), and
    * create BLAST databases for the uPEP databases

Once complete:
    * the new databases will be pushed to production (update the old), and
    * the RefSeq databses counter will be updated 


Usage::

    $ fab -l

    fabfile for doing uPEP things

    Mitchell Stanton-Cook
    m.stantoncook@gmail.com

    TODO: No real verification of the successful completion of step
    (i.e. verify all files downloaded, verify compaction, compilation,
    uPEPFinding & BLAST DB creation). 


    Available commands:

        build_upep_dbs           Update the uPEP databases
        get_NCBI_RefSeq_release  Prints & returns the NCBI RefSeq release number
        get_uPEP_RefSeq_release  Prints & returns the uPEP RefSeq release number


    $ fab -d build_upep_dbs

    Displaying detailed information for task 'build_upep_dbs':

        Update the uPEP databases
        
        Alternatively if a key is given (one of):
            * RefSeq-complete
            * RefSeq-fungi
            * RefSeq-invertebrate
            * RefSeq-plant
            * RefSeq-vertebrate_mammalian
            * RefSeq-vertebrate_other
        
        the task will only upgrade the givenDB.
        
        You can set the outpath with outpath='/dump/me/here'
        
        Setting override to True will update even in local version is the same as
        the NCBI RefSeq version
        
        This task:
            * checks that you need to update (unless overridden)
            * sets up a download directory
            * fetches all required RefSeq databases from NCBI
            * compacts
            * compiles
            * uPEP finds
            * converts to BLAST db
            * moves the files to the production server (creating correct 
              permissions, symbolic links etc)
            * stores the new local RefSeq version
        
        :param outpath: the base location to dump the files to (must exist)
                        i.e. /var/RefSeq/staging (full path as a string). If not 
                        given will be dumped in the fabfile directory
        :param key: [def=None] a specific database to grab. If not given will grab 
                    all
        :param override: [default = False] upgrade even if local and remote are 
                         the same version


    $ fab build_upep_dbs

    Remote Release 59
    Local Release 59
    No updrade required


    $ fab build_upep_dbs:key=RefSeq-plant,override=True

    Remote Release 59
    Local Release 59
    Working with database RefSeq-plant
    Compacting plant.2.rna.gbff.gz
    <SNIP SNIP SNIP>



TODO
----

Future improvements:
    * No real verification of the successful completion of a step. Should
      verify: all files downloaded, compaction, compilation, uPEP finding  
      and BLAST DB creation)
    * When pushing the new uPEP databases to production, the uPEP web app
      should be put into maintainence mode.

