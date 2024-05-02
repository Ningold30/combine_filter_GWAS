import subprocess

#plotting sum stats
def plotting_man_qqplot(input_file
			,input_dir,title):
	print("--Plot set to True")
	plot_script_path= "/mnt/backedup/home/nathanI/Scripts/combine_filter_GWAS/plot/qsub_ManQQplot.sh"
	subprocess.run(['qsub', '-v', f'file1={input_file},dirR={input_dir},title={title}', plot_script_path])
	print("plotting")
	subprocess.run(['qstat','-u nathanI'])
