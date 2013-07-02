RandAL
======

Randomized Sequence Aligner

1. Indexing a reference genome:
randal-build ref_file_name.fasta [index_file_name]

2. Aligning reads to the reference:
randal index_file_name read_file_name.fq map_file_name.sam [OPTIONS]