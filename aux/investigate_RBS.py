#Very basic idea of how to optimize RBS.
#More clever optimization might utilize genetic algorithms.

#>infA_context
#GATAAGGAATTTTTCGCGTCAGGTAACGCCCATCGTTTATCTCACCGCTCCCTTATACGT
#TGCGCTTTTGGTGCGGCTTAGCCGTGTGTTTTCGGAGTAATGTGCCGAACCTGTTTGTTG
#CGATTTAGCGCGCAAATCTTTACTTATTTACAGAACTTCGGCATTATCTTGCCGGTTCAA

#                             *start*
#ATTACGGTAGTGATACCCCAGAGGATTAGATGGCCAAAGAAGACAATATTGAAATGCAAG
#GTACCGTTCTTGAAACGTTGCCTAATACCATGTTCCGCGTAGAGTTAGAAAACGGTCACG
#TGGTTACTGCACACATCTCCGGTAAAATGCGCAAAAACTACATCCGCATCCTGACGGGCG
#ACAAAGTGACTGTTGAACTGACCCCGTACGACCTGAGCAA

import subprocess
import os
import random
import math
import sys

#New users: change these paths.
os.environ['NUPACKHOME'] = "/home/tb/Tools/nupack3.0.5/" #'/home/tb/NuPACK/nupack3.0.4/'
os.environ['PATH'] = os.environ['PATH'] + ':' + "/home/tb/Tools/nupack3.0.5/bin/" #'/home/tb/NuPACK/nupack3.0.4/bin/'
rbs_calc = "/home/tb/Tools/Ribosome-Binding-Site-Calculator-v1.0-master/Run_RBS_Calculator.py"

def score_rbs(seq):
	process = subprocess.Popen(["python", rbs_calc, seq], stdout=subprocess.PIPE)
	out, err = process.communicate()
	return out

def main(thread):
	prefix =  "ctcgaccgtgtcgatgccgat"
	suffix = "gaatggcgaaaga" #sequence is variable beyond this point.
	counter = 0
	out_file = open("rbs_report_" + str(thread) + ".txt", "w", 100)
	#for v0 in ['gtt', 'gtc', 'gta', 'gtg']:
		#for d1 in ['gat', 'gac']:
			#for a0 in ['gct', 'gcc', 'gca', 'gcg']:
				#for d2 in ['gat', 'gac']:
	if True:
		if True:
			if True:
				if True:
					for v in ['gtt', 'gtc', 'gta', 'gtg']:
						for e in ['gaa', 'gag']:
							for y in ['tac', 'tat']:
								for q in ['caa', 'cag']:
									for a1 in ['gct', 'gcc', 'gca', 'gcg']:
										for a2 in ['gct', 'gcc', 'gca', 'gcg']:
											for a3 in ['gct', 'gcc', 'gca', 'gcg']:
												for s in ['agt', 'agc', 'tct', 'tcc', 'tca', 'tcg']:
													counter += 1
													if counter % 8 == thread:
														full = prefix + v + e + y + q + a1 + a2 + a3 + s + suffix
														result = score_rbs(full)
														out_file.write(">" + str(counter) + " " + full + result.split("\n47")[1].split("\n")[0] + "\n")
	out_file.close()

thread = int(sys.argv[-1])
main(thread)
