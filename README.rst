uPEP-helpers
============

A core fabfile for managing uPEP databases

Mitchell Stanton-Cook (m.stantoncook@gmail.com)

Requirements (python libraries):
    * ftputil
    * fabric

Requirements (3rd party):
    * blast2 (legacy)

Overview
--------

Aids in updating uPEP web app databases. If new RefSeq release will:
    * download all required NCBI RefSeq databases, 
    * compact the databases, 
    * compile the databases, 
    * apply uPEP finding (creating uPEP databases), and
    * create BLAST databases for the uPEP databases

Once complete:
    * the new databases will be pushed to production (update the old), and
    * the RefSeq databses counter will be updated 


Usage::

    fab -l
    fab -D build_upep_dbs
    fab build_upep_dbs
    fab build_upep_dbs:key='RefSeq-fungi'


TODO
----

Future improvements:
    * No real verification of the successful completion of a step. Should
      verify: all files downloaded, compaction, compilation, uPEP finding  
      and BLAST DB creation)
    * When pushing the new uPEP databases to production, the uPEP web app
      should be put into maintainence mode.

