from importlib.resources import files
import os
import argparse
import subprocess
import sys
import pandas as pd

def bash_command(cmd):
	p = subprocess.Popen(cmd, shell=True)
	while True:
		return_code = p.poll()
		if return_code is not None:
			break
	return

def run_minimap2(type, reads, ref, prefix, cpu):
	cmd = f"mkdir -p minimap2"
	bash_command(cmd)
	cmd = f"mkdir -p summary"
	bash_command(cmd)
	if type == "paired":
		r1, r2 = reads
		cmd = f"minimap2 -t {cpu} -ax sr {ref} {r1} {r2} > minimap2/{prefix}.bam"
		bash_command(cmd)
	elif type == "single":
		r1 = reads[0]
		cmd = f"minimap2 -t {cpu} -ax sr {ref} {r1} > minimap2/{prefix}.bam"
		bash_command(cmd)
	elif type == "ont":
		r1 = reads[0]
		cmd = f"minimap2 -ax map-ont {ref} {r1} > minimap2/{prefix}.bam"
	cmd = f"samtools sort minimap2/{prefix}.bam > minimap2/{prefix}.sorted.bam"
	bash_command(cmd)
	cmd = f"samtools view -b -F 2308 minimap2/{prefix}.sorted.bam > minimap2/{prefix}.sequins.bam"
	bash_command(cmd)
	cmd = f"samtools index minimap2/{prefix}.sequins.bam"
	bash_command(cmd)
	cmd = f"samtools flagstat minimap2/{prefix}.sorted.bam > summary/{prefix}.mapping.stats"
	bash_command(cmd)
	cmd = f"rm minimap2/{prefix}.bam"
	bash_command(cmd)
	if os.path.exists(f"./minimap2/{prefix}.sequins.bam"):
		return f"minimap2/{prefix}.sequins.bam"
	else:
		return f"error, something wrong with mapping steps"

def isoquant(bam, gtf, ref, cpu, type, prefix):
	output = f"./{prefix}_isoquant"
	cmd = f"isoquant.py --reference ${ref} --genedb ${gtf} --bam {bam} --data_type {type} -o ${output} -t {cpu} --count_exons"
	bash_command(cmd)
	if os.path.exists(f"./{prefix}_isoquant/OUT/OUT.transcript_counts.tsv"):
		return f"./{prefix}_isoquant/OUT"
	else:
		return f"error, isoquant may fail, please check the log in isoquant output folder"

def run_salmon(type, reads, trans_ref, prefix, cpu):
	#create ref index
	cmd = f"salmon index -t ${trans_ref} -i sequins_index -k 31"
	bash_command(cmd)
	if type == "paired":
		r1, r2 = reads
		cmd = f"salmon quant -i sequins_index -l IU -1 {r1} -2 {r2} --validateMappings -p {cpu} -o {prefix}_salmon"
		bash_command(cmd)
	elif type == "single":
		r = reads[0]
		cmd = f"salmon quant -i sequins_index -l SF -r {r} --validateMappings -p {cpu} -o {prefix}_salmon"
		bash_command(cmd)
	else:
		return f"error, unsuppport type"
	if os.path.exists(f"./{prefix}_salmon/quant.sf"):
		return f"./{prefix}_salmon/quant.sf"
	else:
		return f"error, salmon may fail, please check the run folder"
	#salmon quant -i gene_index -l SF -r ${fastq_home}/${b}.fastq --validateMappings -p 20 -o ${b}_quant

def calculate_correlation(quant_out, quant_type, valid_table):
	df_quant = ""
	if quant_type == "isoquant":
		df_count = pd.read_csv(f"{quant_out}/OUT.transcript_counts.tsv", sep="\t")
		df_count.columns = ["Trans_ID", "ISOQUANT_Counts"]
		df_tpm = pd.read_csv(f"{quant_out}/OUT.transcript_tpm.tsv", sep="\t")
		df_tpm.columns = ["Trans_ID", "ISOQUANT_Tpm"]
		df_quant = df_count.merge(df_tpm, on="Trans_ID", how="left")
	elif quant_type == "salmon":
		df_quant = pd.read_csv(quant_out, sep="\t")[["Name","NumReads","TPM"]]
		df_quant.columns = ["Trans_ID", "Salmon_Counts", "Salmon_Tpm"]
	#calculate correlation
	sequin_df = pd.read_csv(valid_table, sep="\t")[["NAME","MIX_A"]]
	sequin_df.columns = ["Trans_ID","Sequin_TRUTH"]
	df_corr = sequin_df.merge(df_quant, on="Trans_ID")
	df_corr = df_corr.drop(columns="Trans_ID")
	corr_result_pearson = df_corr.corr(method="pearson")
	columns = corr_result_pearson.columns
	corr_result_pearson.set_index(columns, inplace=True)
	corr_result_pearson.to_csv(f"./summary/pearson_correlation.csv")
	corr_result_spearman = df_corr.corr(method="spearman")
	corr_result_spearman.set_index(columns, inplace=True)
	corr_result_spearman.to_csv(f"./summary/spearman_correlation.csv")


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--r1",  help="single reads or R1 in paired reads")
	parser.add_argument("--r2",  help="optional, R2 in paired reads")
	parser.add_argument("--ont", action='store_true', help="setting for Long reads Nanopore Sequence")
	parser.add_argument("--work_dir", "-w", default=".", help="work path to store all outputs")
	parser.add_argument("--prefix", default="sequins_qc", help="prefix for all the output")
	parser.add_argument("--cpu", default=5, help="number of threads/cpus to assign to the task")
	args = parser.parse_args()
	#checking sequins data from the resources
	gtf_file = files("sequins_qc").joinpath("data/rnasequin_annotation_2.4.gtf")
	wg_fa = files("sequins_qc").joinpath("data/rnasequin_decoychr_2.4.fa")
	trans_fa = files("sequins_qc").joinpath("data/rnasequin_sequences_2.4.fa")
	valid_table = files("sequins_qc").joinpath("data/sequin_expression.tab")

	#checking tools
	#minimap2
	#isoquant

	#go to the work dir
	try:
		os.chdir(args.work_dir)
	except:
		print(f"no permission in {args.work_dir}, will keep working in your current dir")

	stype = ""
	reads_list = []
	if args.ont:
		stype = "ont"
		reads_list.append(args.r1)
	elif args.r2:
		stype = "paired"
		reads_list.append(args.r1)
		reads_list.append(args.r2)
	else:
		stype = "single"
		reads_list.append(args.r1)
	#mapping
	path_to_bam = run_minimap2(stype, reads_list, wg_fa, args.prefix, args.cpu)
	if "error" in path_to_bam:
		print(path_to_bam)
		sys.exit()
	quant_out = ""
	quant_type = ""
	if args.ont:
		#quant
		quant_type = "isoquant"
		quant_out = isoquant(path_to_bam, gtf_file, wg_fa, args.cpu, "nanopore",args, prefix)
	else:
		quant_type = "salmon"
		quant_out = run_salmon(stype, reads_list, trans_fa, args.prefix, args.cpu)

	if not "error" in quant_out:
		calculate_correlation(quant_out, quant_type, valid_table)
		print("Wonderful, sequins_qc finished, please check summary table in summary folder")
	else:
		print("Error in quantification steps, please check output file to debug")
