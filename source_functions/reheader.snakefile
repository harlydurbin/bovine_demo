# snakemake -s source_functions/reheader.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200706.reheader.log

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
	 	expand("data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz", chr = config['chr'])

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
