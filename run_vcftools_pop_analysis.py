#!/bin/env python3

import os, sys, argparse, subprocess

def main():

	## arguments with help messages
	parser = argparse.ArgumentParser(description="Running vcftools Fst windowed, diversity, snp_density. Need 2 pop files - can be females and males.", allow_abbrev=False)
	parser.add_argument("--vcf_files", nargs = "*", help = "vcf file(s) to use. Default: No", default = []) #default = empty list
	parser.add_argument("--pop1_file", type = str, default = "",help = "File containing the individuals for population 1. Default: No")
	parser.add_argument("--pop2_file", type = str, default = "", help = "File containing the individuals for population 2. Default: No")
	parser.add_argument("--out_dir", type = str, default = "", help = "Ouptut directory for the pop stat files. Default: No")
	parser.add_argument("--prefix_out", type = str, help = "Output prefix. Default: No")
	parser.add_argument("--wind_size", type = int, help = "Window size for the windowed analysis. Default: 10000", default = 10000)

	## Parsing the arguments
	args = parser.parse_args()
#	print("Arguments: ",args,"\n")
	
	for filename in args.vcf_files:
		if filename.endswith(".vcf.gz"):
			file_base=os.path.basename(filename)
			file_no_gz = os.path.splitext(file_base)[0]
			file_no_ext = os.path.splitext(file_no_gz)[0]
			file_base_pop1 = os.path.basename(args.pop1_file)
			file_no_ext_pop1 = os.path.splitext(file_base_pop1)[0]
			file_base_pop2 = os.path.basename(args.pop2_file)
			file_no_ext_pop2 = os.path.splitext(file_base_pop2)[0]
			cmd = f"vcftools --gzvcf {filename} --weir-fst-pop {args.pop1_file} --weir-fst-pop {args.pop2_file} --out {args.out_dir}/{file_no_ext}_{file_no_ext_pop1}_VS_{file_no_ext_pop2} --fst-window-size {args.wind_size}"
			print(f"Command: {cmd}")
			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
			process.wait()
# 			cmd = f"vcftools --gzvcf {filename} --keep {args.pop1_file} --out {args.out_dir}/{file_no_ext}_{file_no_ext_pop1} --window-pi {args.wind_size}"
# 			print(f"Command: {cmd}")
# 			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
# 			process.wait()
# 			cmd = f"vcftools --gzvcf {filename} --keep {args.pop2_file} --out {args.out_dir}/{file_no_ext}_{file_no_ext_pop2} --window-pi {args.wind_size}"
# 			print(f"Command: {cmd}")
# 			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
# 			process.wait()
#			cmd = f"vcftools --gzvcf {filename} --keep {args.pop1_file} --out {args.out_dir}/{file_no_ext}_{file_no_ext_pop1} --SNPdensity {args.wind_size}"
#			print(f"Command: {cmd}")
#			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
#			process.wait()
#			cmd = f"vcftools --gzvcf {filename} --keep {args.pop2_file} --out {args.out_dir}/{file_no_ext}_{file_no_ext_pop2} --SNPdensity {args.wind_size}"
#			print(f"Command: {cmd}")
#			process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines = True, shell = True)
#			process.wait()			
if __name__ == "__main__":
	main()
