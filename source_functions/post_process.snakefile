# snakemake -s source_functions/post_process.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200728.joint_genotyping.log

# include path is relative to the path of this file
include: "joint_genotyping.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

import os

os.makedirs("log/slurm_out/post_process", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/post_process/{rules}", rules = config['post_process_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['post_process_rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("temp/post_process", exist_ok = True)
for x in expand("temp/post_process/{rules}", rules = config['post_process_rules']):
    os.makedirs(x, exist_ok = True)

rule post_process_all:
	input:
	 	expand("data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz", chr = config['chr']), expand("data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz.tbi", chr = config['chr']), expand("data/derived_data/joint_genotyping/collect_metrics/collect_metrics.{chr}.variant_calling_detail_metrics", chr = config['chr']), expand("data/derived_data/joint_genotyping/collect_metrics/collect_metrics.{chr}.variant_calling_summary_metrics", chr = config['chr']),
		expand("data/derived_data/joint_genotyping/validate_variants/validate_variants.{chr}.txt", chr = config['chr'])

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
		java_tmp = "temp/post_process/remove_samples/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/remove_samples/remove_samples.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -L {params.chr} -T SelectVariants -nt {params.nt} -R {params.ref_genome} -xl_sf {input.remove_list} -ef --removeUnusedAlternates -env -V {input.vcf} -o {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule collect_metrics:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz.tbi"
	params:
		picard_module = config['picard_module'],
		java_tmp = "temp/joint_genotyping/collect_metrics/{chr}",
		gc_threads = config['collect_metrics_gc'],
		xmx = config['collect_metrics_xmx'],
		picard_path = config['picard_path'],
		truth_file = config['truth_file'],
		ref_genome = config['ref_genome'],
		nt = config['collect_metrics_nt'],
		psrecord = "log/psrecord/joint_genotyping/collect_metrics/collect_metrics.{chr}.log",
		prefix = "data/derived_data/joint_genotyping/collect_metrics/collect_metrics.{chr}"
	output:
		detail = "data/derived_data/joint_genotyping/collect_metrics/collect_metrics.{chr}.variant_calling_detail_metrics",
		summary = "data/derived_data/joint_genotyping/collect_metrics/collect_metrics.{chr}.variant_calling_summary_metrics"
	shell:
		"""
		module load {params.picard_module}
		psrecord "java -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.picard_path} CollectVariantCallingMetrics TMP_DIR={params.java_tmp} INPUT={input.vcf} OUTPUT={params.prefix} DBSNP={params.truth_file} R={params.ref_genome} THREAD_COUNT={params.nt}" --log {params.psrecord} --include-children --interval 5
		"""

rule validate_variants:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{chr}.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['validate_variants_gc'],
		xmx = config['validate_variants_xmx'],
		java_tmp = "temp/joint_genotyping/validate_variants/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/validate_variants/validate_variants.{chr}.log",
		chr = "{chr}"
	output:
		report = "data/derived_data/joint_genotyping/validate_variants/validate_variants.{chr}.txt"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T ValidateVariants -R {params.ref_genome} -V {input.vcf} -L {params.chr} --warnOnErrors &> {output.report}" --log {params.psrecord} --include-children --interval 5
		"""
