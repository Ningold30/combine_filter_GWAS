#!/usr/bin/env python

##########################
#import
##########################
import pandas as pd
import os
import argparse
import numpy as np
import subprocess
from munge.munging import munge_results
from clump.clumping import clumping_results
from plot.QQMan_PLotting import plotting_man_qqplot
import pathlib

##########################
#define arguments
##########################

#combining, filtering  and processing sumstats
def combine_and_filter(file_chr_22
			, output_file
			, regenie=False
			, saige_binary=False
			, saige_quant=False
			, min_af=0.01
			, min_info=0.3
			, num_chrom=22
			, clump=False
			, clumpr2=0.001
			, clumpKB=1000
			, clumpP=5e-8
			, munge=False
			, ManQQplot=False
			, ManQQplot_title="Title"
			):
	print(f"using follwoing parameters:\nfile_chr_22 ={file_chr_22}\noutput_file={output_file}\nmin_af={min_af}\nnum_chrom={num_chrom}\nmunge={munge}\nManQQplot={ManQQplot}\nManQQplot_title={ManQQplot_title}\nclump={clump}\nclumpr2={clumpr2}\nclumpKB={clumpKB}\nclumpP={clumpP}")
	# Create empty DataFrame
	combined_df = pd.DataFrame()
	# deconsruct the file name
	input_folder,file_chr_22_name = os.path.split(file_chr_22)
	prefix,suffix = file_chr_22_name.split("22")
	# Loop through the specified number of chromosomes
	for chrom in range(1, num_chrom + 1):
		print("for chromosome",chrom)
		# Construct the filename using pre and suffix
		filename = os.path.join(input_folder, f"{prefix}{chrom}{suffix}")
		print(f"Reading file: {filename}")
		# Read the GWAS summary statistics file for each chromosome
		df = pd.read_csv(filename, sep = "\s+")
		print("Columns of file\n",df.columns)
		#make columns the same accross softwares
		if regenie:
			print("summstats format set using --regenie")  
			#change names A1 is the effect allele, ALLELE1 in regenie
			df.columns = ["CHR", "BP", "SNP", "A2", "A1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P", "EXTRA"]
			# convert log10 pval to pvalue
			df['P'] = 10 ** (-df['LOG10P'])
			df = df.loc[:, ["CHR", "BP", "SNP", "A1", "A2", "A1FREQ", "INFO", "N", "BETA", "SE", "P"]]
			print("new columns",df.columns)
		
		if saige_binary:
			print("summstats format set using --saige_binary") 
			#change column names, A1 is the effect alleel column (Allele2 in saige)
			df.columns = ["CHR", "BP", "SNP", "A2", "A1", "AC_Allele2", "A1FREQ", "MissingRate","BETA", "SE", "Tstatvar", "P", "p.value.NA", "Is.SPA", "AF_case", "AF_ctrl", "N_case", "N_ctrl", "N_case_hom", "N_case_het", "N_ctrl_hom", "N_ctrl_het"]
			#get analagous to INFO score
			df['INFO'] = 1-df['MissingRate']
			#get N column
			df['N'] = df['N_case'] + df['N_ctrl']
			df = df.loc[:, ["CHR", "BP", "SNP", "A1", "A2", "A1FREQ", "INFO", "N", "BETA", "SE", "P"]]
		if saige_quant:
			print("summstats format set using --saige_quant")
			#change column names, A1 is the effect alleel column (Allele2 in saige)
			df.columns = ["CHR", "BP", "SNP", "A2", "A1", "AC_Allele2", "A1FREQ", "MissingRate","BETA", "SE", "Tstat", "var", "P", "N"]
			df['INFO'] = 1-df['MissingRate']
			df = df.loc[:, ["CHR", "BP", "SNP", "A1", "A2", "A1FREQ", "INFO", "N", "BETA", "SE", "P"]]
		# Filter rows based on AF and INFO
		print("Line count before filtering",df.shape[0])
		df = df[(df['A1FREQ'] >= min_af) & (df['A1FREQ'] <= 1-min_af) & (df['INFO'] >= min_info)]
		print(f"Line count after filtering:, AF {min_af}, INFO {min_info}", df.shape[0])
		#ammend current df to combined_df
		combined_df = pd.concat([combined_df, df], ignore_index=True)
		print("Line count of combined:", combined_df.shape[0])
		print(combined_df.head)
	#save
	output_file = f"{output_file}_INFO_{min_info}_MAF_{min_af}.txt"
	output_name= os.path.join(input_folder, output_file)
	combined_df.to_csv(output_name, index=False, sep = " ", na_rep ="NA")
	print(f"File saved to: {output_name}")
	
	if munge:
		#format for LDSC
		#make LDSC folder
		LDSC_output_dir = f"{input_folder}/LDSC"
		pathlib.Path(LDSC_output_dir).mkdir(parents=True, exist_ok=True)
		#os.mkdir(LDSC_output_dir, exist_ok=True)
		#subset LDSC columsn
		LDSC = combined_df.loc[:, ["SNP","A1","A2","BETA","P","N"]]
		#save
		LDSC_output = f"{LDSC_output_dir}/{output_file}.LDSC"		
		LDSC.to_csv(f"{LDSC_output}",index=False, sep = " ", na_rep = "NA")
		print(f"LDSC files saved here: {LDSC_output}")
		#LDSC_path, LDSC_file = os.path.split(LDSC_output)
		munge_results(f"{output_file}.LDSC", LDSC_output_dir)
		
	if ManQQplot:
		plotting_man_qqplot(output_file,input_folder,ManQQplot_title)

	if clump:
		CLUMP = combined_df.loc[:, ["SNP","P"]]
		print("--clump set to true")
		print(f"clumping with the following thresholds: R2={clumpr2}, P value={clumpP}, Kilabases={clumpKB}")
		temp_clump = f"{output_file}.CLUMP.TMP"
		CLUMP.to_csv(f"{input_folder}/{temp_clump}",index=False, sep = " ", na_rep = "NA")
		clumping_results(input_file = temp_clump
				,input_dir = input_folder
				,pval = "P"
				,out_file = output_file
				,r2 = clumpr2
				,KB = clumpKB
				,P = clumpP
				)		

###########################
#Parse arguments	
###########################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="This Script combines and filters GWAS summary statistics across chromsomes, allows you to munge the summary statistics (--munge), allows you to plot a Mannhattan plot and qqplot using the --ManQQplot and --ManQQplot_title, and clump (--clump, --clumpr2, --clumpP, --clumpKB")
	parser.add_argument("--file_chr_22", required=True, help="the full path and  file name of the results of chromsome 22 - used to work out naming convention of input and output files")
	parser.add_argument("--output_file", required=True, help="name for output file")
	parser.add_argument("--regenie", required=False, action="store_true", default=False, help="was regenie used to generate these results")
	parser.add_argument("--saige_binary", required=False, action="store_true", default=False, help="was saige used to perform this GWAS as a binary trait")
	parser.add_argument("--saige_quant", required=False, action="store_true", default=False, help="was saige used to perform this GWAS as a quantitiative trait")
	parser.add_argument("--min_af", required=False, type=float, default=0.01, help="minimum allele frequency, default = 0.01")
	parser.add_argument("--min_info", required=False, type=float, default=0.3, help="minimum info score, default = 0.3")
	parser.add_argument("--num_chrom", required=False, type=int, default=22, help="number of chromosomes, default = 22" )
	parser.add_argument("--clump", required=False, action="store_true", default=False, help="Do you want to clump the results?")
	parser.add_argument("--clumpr2", required=False, type=float, default=0.001, help="r2 threshold used for clumping, default 0.001")
	parser.add_argument("--clumpKB", required=False, type=int, default=1000, help="kilobases used for clumping, default 1000KB")
	parser.add_argument("--clumpP", required=False, type=float, default=5e-8, help="genome-wide significant P value, used for clumping, default 5e-8")
	parser.add_argument("--munge", required=False, action="store_true", default=False, help="do you want to munge the sumstats for LDSC, default = False")
	parser.add_argument("--ManQQplot", required=False, action="store_true", default=False, help="do you want to plot a Manhattan and QQ plt, default = False, must be accompanied with the --ManQQplot_title flag")
	parser.add_argument("--ManQQplot_title",type=str, required=False, default="title", help="provide a title for the plots") 

	args = parser.parse_args()
	
	# Call the function to combine and filter the data
	combine_and_filter(args.file_chr_22
			, args.output_file
			, args.regenie
			, args.saige_binary
			, args.saige_quant
			, args.min_af
			, args.min_info
			, args.num_chrom
			, args.clump
			, args.clumpr2
			, args.clumpKB
			, args.clumpP
			, args.munge
			, args.ManQQplot
			, args.ManQQplot_title
			)






