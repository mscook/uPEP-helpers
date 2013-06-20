uPEP-helpers
============

A core fabfile for managing uPEP databases

Mitchell Stanton-Cook (m.stantoncook@gmail.com)

Requirements (python libraries):
    * ftputil
    * fabric

Requirements (3rd party):
    * blast (legacy)

Usage::

    fab -l
    fab -D build_upep_dbs
    fab build_upep_dbs
    fab build_upep_dbs:key='RefSeq-fungi'

