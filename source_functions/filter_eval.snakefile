# snakemake -s source_functions/filter_eval.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --resources load=100 --config --cluster-config source_functions/cluster/filter_eval.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200720.filter_eval.log

import os

configfile: "source_functions/config/filter_eval.config.yaml"

os.makedirs("log/slurm_out/filter_eval", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/filter_eval/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("temp/filter_eval/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
	 	expand("data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.table", dataset = config['dataset']), expand("data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.ldepth.mean", dataset = config['dataset']), "data/derived_data/joint_genotyping/filter_eval/filter_eval.post.table", expand("data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.GQ.FORMAT", dataset = config['dataset'])



def input_chooser(WC):
	# References the dictionary, then the value associated with that key (which is your wildcard)
	file_location = config['dataset_input'][WC.dataset]
	return file_location

# Create a chr28 vcf with snps & indels to evaluate how many indels are within
# 5bp of an inel, filtering values
rule biallelic_28:
	input:
		vcf = input_chooser
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		nt = config['biallelic_28_nt'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/filter_eval/biallelic_28"
	output:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz.tbi"
	shell:
		"""
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 -jar {params.gatk_path} -nt {params.nt} -T SelectVariants -R {params.ref_genome} -L 28 -V {input.vcf} --restrictAllelesTo BIALLELIC -selectType SNP -o {output.vcf}
		"""

# vcftools --gzvcf data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz --get-INFO GQ --out data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28

rule depth_28:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28"
	output:
		depth = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.ldepth.mean"
	shell:
		"""
		module load vcftools
		vcftools --gzvcf {input.vcf} --site-mean-depth --out {params.prefix}
		"""

rule table_28:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/joint_genotyping/table_28"
	output:
		table = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.table"
	shell:
		"""
		module load {params.java_module}
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 -jar {params.gatk_path} -R {params.ref_genome} -L 28 -T VariantsToTable -V {input.vcf} -F POS -F TYPE -F TRANSITION -F QD -F FS -F MQ -F ReadPosRankSum -F MQRankSum -F NO-CALL -F N-CALLED -F VAR -o {output.table}
		"""

rule table_28_post:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/remove_samples.28.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/joint_genotyping/table_28_post"
	output:
		table = "data/derived_data/joint_genotyping/filter_eval/filter_eval.post.table"
	shell:
		"""
		module load {params.java_module}
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 -jar {params.gatk_path} -R {params.ref_genome} -L 28 -T VariantsToTable -V {input.vcf} -F POS -F TYPE -F TRANSITION -F QD -F FS -F MQ -F ReadPosRankSum -F MQRankSum -F NO-CALL -F N-CALLED -F VAR -o {output.table}
		"""

rule select_random:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/joint_genotyping/select_random",
		fraction = config['fraction'],
		nt = config['biallelic_28_nt']
	output:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 -jar {params.gatk_path} -nt {params.nt} -T SelectVariants -R {params.ref_genome} -L 28 -V {input.vcf} --select_random_fraction {params.fraction} -o {output.vcf}
		"""

rule extract_gq:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random"
	output:
		gq = "data/derived_data/joint_genotyping/filter_eval/filter_eval.{dataset}.28_random.GQ.FORMAT"
	shell:
		"""
		module load vcftools
		vcftools --gzvcf {input.vcf} --extract-FORMAT-info GQ --out {params.prefix}
		"""
