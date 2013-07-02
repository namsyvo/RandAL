#GA_v2.5


0) Change scripts to run experiments "one time"
    - Use data from a fixed path (don't need to copy to current directory).
	- run alignment and evaluation in one script.

1) do_exps.py:

dir_names = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep', 'Soran']
rlen = [600, 800]

with parameters:
rlen = [35, 51, 76, 100, 200, 400]
	dir_name = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep']
		th = int(math.floor(error_prob * length + th_para * math.sqrt(length * error_prob * (1 - error_prob))))
		att = att_para * (th+1)
		minws = int(math.ceil(length/(th + 1))) + 3

		(error_prob, att_para, th_para) = (0.02, 1, 4)

The only thing that is different from version 2.3 is minws = int(math.ceil(length/(th + 1))) + 5

max of candidates is 2
positions to search for randomized: (10, length-10)

Output files for each species and each length:
- *.sam (for all reads and for false reads)
- *.all_time.txt

To evaluate: eval_exps.py


1) do_exps_soran.py: run for Soran

dir_names = ['Soran']
rlen = [35, 51, 76, 100, 200, 400]

with parameters:
rlen = [35, 51, 76, 100, 200, 400]
	dir_name = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep']
		th = int(math.floor(error_prob * length + th_para * math.sqrt(length * error_prob * (1 - error_prob))))
		att = att_para * (th+1)
		minws = int(math.ceil(length/(th + 1))) + 3

		(error_prob, att_para, th_para) = (0.02, 1, 4)

The only thing that is different from version 2.3 is minws = int(math.ceil(length/(th + 1))) + 5

max of candidates is 2
positions to search for randomized: (10, length-10)

Output files for each species and each length:
- *.sam (for all reads and for false reads)
- *.all_time.txt

Take all results to file overall_result.txt


2) do_exps_1.py and do_exps_4.py
Run experiments for 1% and 4% base eroor rates of simulated reads

dir_names = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep', 'Soran']
rlen = [35, 400]

with parameters:
rlen = [35, 51, 76, 100, 200, 400]
	dir_name = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep']
		th = int(math.floor(error_prob * length + th_para * math.sqrt(length * error_prob * (1 - error_prob))))
		att = att_para * (th+1)
		minws = int(math.ceil(length/(th + 1))) + 3

		(error_prob, att_para, th_para) = (0.01, 1, 4) (for 1%)
		(error_prob, att_para, th_para) = (0.04, 1, 4) (for 4%)

read files:
	Dros-35-1.fq and others
	Dros-35-4.fq and others

The only thing that is different from version 2.3 is minws = int(math.ceil(length/(th + 1))) + 5

max of candidates is 2
positions to search for randomized: (10, length-10)

Output files for each species and each length:
- *.sam (for all reads and for false reads)
- *.all_time.txt

Evaluate results:
*accuracy.txt
*running_time.txt

Take all results to file overall_result-1.txt
Take all results to file overall_result-4.txt



2) Evaluation program: wgsim_eval_tofile.pl produces:

File: accuracy.txt
error	aligned reads
each next line for each set of parameters

File: running_time.txt
Time
each next line for each set of parameters


3) Take overrall results: process_results.py
dir_names = ['Dros', 'Ecoli', 'Pseu', 'Stap', 'Strep', 'Soran']
rlen = [600, 800]

output file: overall_result.txt

Species
randomized_w5	error	mapped reads	time	error	mapped reads	time	...
Species
randomized_w5	error	mapped reads	time	error	mapped reads	time	...
...
