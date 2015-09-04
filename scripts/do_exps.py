'''
Run the alignment for species and lengths
Evaluate the alignment results
'''

import os
import sys
import re
import math

def cal_paras(length, error_prob, att_para, th_para):
	th = int(math.floor(error_prob * length + th_para * math.sqrt(length * error_prob * (1 - error_prob))))
	att = att_para * (th + 1)
	#test with ws
	minws = int(math.ceil(length/(th + 1))) + 3
	return th, att, minws

dir_names = ['Drosophila3R.fasta', 'Staphylococcus.fasta']
rlen = [35, 51, 76, 100, 200, 400]

for d in dir_names:
	if not os.path.isdir(d):
		os.mkdir(d)
	for length in rlen:
		if not os.path.isdir(d + '/' + str(length)):
			os.mkdir(d + '/' + str(length))

		th, att, minws = cal_paras(length, 0.02, 1, 4)
		B_file = '../data/genomes/' + d + '/B_100_index'
		F_file = '../data/genomes/' + d + '/F_100_index'
		ref_file = '../data/genomes/' + d + '/' + d
		read_file = '../data/reads/' + d + '/' + d + '-' + str(length) + '.fq'
		result_path = d + '/' + str(length) +'/results/'
		if not os.path.isdir(result_path):
			os.mkdir(result_path)

		eval_path = result_path + 'Eval/'
		if not os.path.isdir(eval_path):
			os.mkdir(eval_path)
		sam_file = result_path + d + '-' + str(length) + '-' + str(th) + '_' + str(att) + '_' + str(minws) + '.sam'
		time_file = eval_path + 'all_time.txt'
		exp_cmd = '(time ./randal-align ' + B_file + ' ' + F_file + ' ' + ref_file + ' ' + read_file + ' ' + sam_file + \
						 ' -th ' + str(th) + ' -att ' + str(att) + ' -mode 3 -rep 10 -minws ' + str(minws) + \
						 ' -maxcdt 2)' + ' 2>>' + time_file + ' &'
		os.system(exp_cmd)
