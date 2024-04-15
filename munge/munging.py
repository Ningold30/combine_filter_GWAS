import subprocess

#munging summstats
def munge_results(input_file,input_dir):
	print("--Munge set to True")
	munge_script_path= "/mnt/backedup/home/nathanI/Scripts/LDSC/munge.sh"
	subprocess.run(['qsub', '-v', f'file1={input_file},dirR={input_dir}', munge_script_path])
	subprocess.run(['qstat','-u nathanI'])
