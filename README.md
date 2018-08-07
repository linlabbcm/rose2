======
ROSE 2
======

Rank ordering of super-enhancers.

Install
=======

```
pip install -e git+https://github.com/linlabbcm/rose2.git
```

Dependencies
============

The CRC software uses the following dependencies:

- ``Bamliquidator``

- ``Samtools``

Usage
=====

As a command line tool:
```
rose2 -g [GENOME] -i [INPUT_REGION_GFF] -r [RANKBY_BAM_FILE] -o [OUTPUT_FOLDER] [OPTIONAL_FLAGS]
```

As a python library:
```
import rose2

rose2.rose(input_file, rankby, output_folder, genome, bams=None, control='', stitch=None, tss=0, mask_file=None)
```
