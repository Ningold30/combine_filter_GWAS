import subprocess

#clumping summ stats
def clumping_results(input_file
			,input_dir
			,pval
			,out_file
			,r2
			,KB
			,P
			):
	clump_script_path="/mnt/backedup/home/nathanI/Scripts/clumping_GWAS/pipeline_clump.sh"
	print(f"function r2={r2}, P value = {P}, KB = {KB}")
	subprocess.run(['qsub', '-v', f'r2={r2},KB={KB},P={P},file1={input_file},dir={input_dir},pval={pval},file2={out_file}', clump_script_path])
	subprocess.run(['qstat','-u nathanI'])
