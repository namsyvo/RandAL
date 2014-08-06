RandAL - Randomized Sequence Aligner
------------------------------------
RandAL is a tool for aligning DNA sequences to reference genomes. The tool is developed based on a new randomized algorithm with the distinction of having high performance across a wide range of read lengths and base error rates.

RandAL is implemented in C++; FM-index codes are adapted from an external library (http://code.google.com/p/fmindex-plus-plus). The tool has been tested with Debian Squeeze 6.0.6 and Mac OS 10.7.5.

More detail about RandAL can be found in our paper "Nam S. Vo, Quang Tran, Nobal Nilaura, Vinhthuy Phan. RandAL: a randomized approach to aligning DNA sequences to reference genomes. BMC Genomics 2014, 15(Suppl 5):S2 doi:10.1186/1471-2164-15-S5-S2". (http://www.biomedcentral.com/1471-2164/15/S5/S2)


Directory organization:
-----------------------

./src: source code of RandAL.
./data: several datasets for testing.
./scripts: scripts to support running and testing RandAL.


RandAL source code:
-------------------

./src
See MANUAL for detailed information on how to use RandAL.
See also LICENSE, VERSION, and CHANGELOG for other information.


Sample datasets for testing RandAL:
-----------------------------------

./data/genomes: store two reference genomes.
	/Staphylococcus.fasta (bacterium, http://www.ebi.ac.uk/ena/data/view/Taxon:663951)
	/Drosophila melanogaster chromosome 3R.fasta (eukaryote, http://www.ebi.ac.uk/ena/data/view/AE014297)
Genomes are taken from EBI (http://www.ebi.ac.uk). Users also can find reference genomes at NCBI (http://www.ncbi.nlm.nih.gov).

./data/reads: store simulated reads for above genomes.
	/Staphylococcus: simulated reads for Staphylococcus.
	/Drosophila3R: simulated reads for Drosophila melanogaster chromosome 3R.
Reads are generated with a simulator named wgsim (https://github.com/lh3/wgsim). 100,000 reads with length from 35bps to 400bps are generated with default settings. See https://github.com/lh3/wgsim for detailed information on how to generate simulated reads and evaluate the alignment results.


Supporting scripts to run and evaluate RandAL:
----------------------------------------------

./scripts/do_exps.py
Python script to test RandAL with multiple simulated datasets:

Usage:  python do_exps.py -r ref -l read_len -e error_rate -o result_file_name
Example: python do_exps.py -r Stap -r Dros -l 35 -l 51 -e 0.02 -e 0.04 -o overall_results.txt

./scripts/wgsim_eval.pl
Perl script to evaluate a SAM output and then produce result to screen. Original code from https://github.com/lh3/wgsim.

./scripts/wgsim_eval_tofile.pl
Perl script to evaluate a SAM ouput and then write results (including error mapped reads) to files, used in the script do_exps.py. Modified from wgsim_eval.pl.


Contact:
--------

nsvo1@memphis.edu
qmtran@memphis.edu
vphan@memphis.edu
