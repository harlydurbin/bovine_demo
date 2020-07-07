# snakemake -s source_functions/reheader.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200707.reheader.log

# include path is relative to the path of this file
include: "joint_genotyping.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

import os

os.makedirs("log/slurm_out/reheader", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/reheader/{rules}", rules = config['reheader_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['reheader_rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("temp/reheader", exist_ok = True)
for x in expand("temp/reheader/{rules}", rules = config['reheader_rules']):
    os.makedirs(x, exist_ok = True)


rule reheader_all:
	input:
	 	expand("data/derived_data/joint_genotyping/reheader/bovine_demo.{chr}.vcf.gz.tbi", chr = config['chr']), expand("data/derived_data/joint_genotyping/reheader/bovine_demo.{chr}.vcf.gz", chr = config['chr'])

rule remove_samples:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz.tbi",
		remove_list = "data/derived_data/joint_genotyping/remove_samples/remove.txt"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['select_variants_gc'],
		nt = config['select_variants_nt'],
		chr = "{chr}",
		java_tmp = "temp/reheader/remove_samples/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/remove_samples/remove_samples.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -L {params.chr} -T SelectVariants -nt {params.nt} -R {params.ref_genome} -xl_sf {input.remove_list} -V {input.vcf} -o {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule sample_list:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.28.vcf.gz"
	params:
		bcftools_module = config['bcftools_module']
	output:
		sample_list = "data/derived_data/joint_genotyping/reheader/sample_list.txt"
	shell:
		"""
		module load {params.bcftools_module}
		bcftools query -l {input.vcf} > {output.sample_list}
		"""

rule sample_key:
	input:
		sample_list = "data/derived_data/joint_genotyping/reheader/sample_list.txt",
		script = "source_functions/sample_key.R"
	params:
		r_module = config['r_module']
	output:
		sample_key = "data/derived_data/joint_genotyping/reheader/sample_key.txt"
	shell:
		"""
		module load {params.r_module}
		Rscript --vanilla {input.script}
		"""

rule reheader:
	input:
		sample_key = "data/derived_data/joint_genotyping/reheader/sample_key.txt",
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz.tbi"
	params:
		bcftools_module = config['bcftools_module'],
		nt = config['reheader_nt'],
		psrecord = "log/psrecord/joint_genotyping/reheader/reheader.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/reheader/bovine_demo.{chr}.vcf.gz"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools reheader --samples {input.sample_key} -O z -o {output.vcf} {input.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule index_reheader:
	input:
		vcf = "data/derived_data/joint_genotyping/reheader/bovine_demo.{chr}.vcf.gz"
	params:
		bcftools_module = config['bcftools_module']
	output:
		tbi = "data/derived_data/joint_genotyping/reheader/bovine_demo.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.bcftools_module}
		tabix {input.vcf}
		"""
