#!/bin/bash
PYTHONPATH=/hpc/local/CentOS7/cog/lib/python2.7/site-packages/:$PYTHONPATH

/hpc/local/CentOS7/common/lang/python/2.7.10/bin/python scripts/sam.py $@

