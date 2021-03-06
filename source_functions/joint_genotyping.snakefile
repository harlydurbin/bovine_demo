# snakemake -s source_functions/joint_genotyping.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200703.joint_genotyping.log

# paste(c(1:29, "X", "Y"), collapse = "', '")

# paste(c(1:139), collapse = "', '")

import os

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/joint_genotyping", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/joint_genotyping/{rules}", rules = config['joint_genotyping_rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/joint_genotyping", exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['joint_genotyping_rules']):
    os.makedirs(x, exist_ok = True)

# Make temp directories if they don't exist
for x in expand("temp/joint_genotyping/{rules}", rules = config['joint_genotyping_rules']):
    os.makedirs(x, exist_ok = True)

rule joint_genotyping_all:
	input:
		expand("data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz", chr = config['chr']), expand("data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz.tbi", chr = config['chr'])

rule combine_gvcfs:
# Can't parallelize CombineGVCFs!!!
	input:
		list = "data/derived_data/joint_genotyping/genotyping_cohorts/{chr}/cohort_{cohort}.list"
	params:
		java_tmp = "temp/joint_genotyping/combine_gvcfs/{chr}/cohort_{cohort}/",
		cohort = "{cohort}",
		chr = "{chr}",
		ref_genome = config['ref_genome'],
		gatk_path = config['gatk_path'],
		java_module = config['java_module'],
		gc_threads = config['combine_gvcfs_gc'],
		xmx = config['combine_gvcfs_xmx'],
		psrecord = "log/psrecord/joint_genotyping/combine_gvcfs/combine_gvcfs.cohort_{cohort}.{chr}.log"
	output:
		gvcf = temp("data/derived_data/joint_genotyping/combine_gvcfs/{chr}/combine_gvcfs.cohort_{cohort}.{chr}.g.vcf.gz"),
		tbi = temp("data/derived_data/joint_genotyping/combine_gvcfs/{chr}/combine_gvcfs.cohort_{cohort}.{chr}.g.vcf.gz.tbi")
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -T CombineGVCFs -R {params.ref_genome} -L {params.chr} -V {input.list} -o {output.gvcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule genotype_gvcfs_list:
	input:
		gvcfs = expand("data/derived_data/joint_genotyping/combine_gvcfs/{{chr}}/combine_gvcfs.cohort_{cohort}.{{chr}}.g.vcf.gz",
		cohort = config['cohort']),
		tbis = expand("data/derived_data/joint_genotyping/combine_gvcfs/{{chr}}/combine_gvcfs.cohort_{cohort}.{{chr}}.g.vcf.gz.tbi",
		cohort = config['cohort'])
	params:
		chr = "{chr}"
	output:
		list = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.list"
	shell:
		"ls -d data/derived_data/joint_genotyping/combine_gvcfs/{params.chr}/combine_gvcfs.*.g.vcf.gz > {output.list}"

rule genotype_gvcfs:
	input:
		list = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.list",
		gvcfs = expand("data/derived_data/joint_genotyping/combine_gvcfs/{{chr}}/combine_gvcfs.cohort_{cohort}.{{chr}}.g.vcf.gz",
		cohort = config['cohort']),
		tbis = expand("data/derived_data/joint_genotyping/combine_gvcfs/{{chr}}/combine_gvcfs.cohort_{cohort}.{{chr}}.g.vcf.gz.tbi",
		cohort = config['cohort'])
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['genotype_gvcfs_gc'],
		xmx = config['genotype_gvcfs_xmx'],
		nt = config['genotype_gvcfs_nt'],
		chr = "{chr}",
		java_tmp = "temp/joint_genotyping/genotype_gvcfs/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -nt {params.nt} -T GenotypeGVCFs -R {params.ref_genome} -L {params.chr} -V {input.list} --useNewAFCalculator --heterozygosity 0.0033 --standard_min_confidence_threshold_for_calling 10 -o {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

# Restrict to biallelic SNPs
# Remove indels
rule select_variants:
	input:
		vcf = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.{chr}.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['select_variants_gc'],
		# xmx = config['select_variants_xmx'],
		nt = config['select_variants_nt'],
		chr = "{chr}",
		java_tmp = "temp/joint_genotyping/select_variants/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/select_variants/select_variants.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/select_variants/select_variants.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/select_variants/select_variants.{chr}.vcf.gz.tbi"
	# -Xmx{params.xmx}g
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -T SelectVariants -nt {params.nt} -R {params.ref_genome} -L {params.chr} -selectType SNP --restrictAllelesTo BIALLELIC -V {input.vcf} -o {output.vcf}" --log {params.psrecord} --include-children --interval 2
		"""

rule variant_filtration:
	input:
		vcf = "data/derived_data/joint_genotyping/select_variants/select_variants.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/select_variants/select_variants.{chr}.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['variant_filtration_gc'],
		xmx = config['variant_filtration_xmx'],
		chr = "{chr}",
		java_tmp = "temp/joint_genotyping/variant_filtration/{chr}",
		gatk_path = config['gatk_path'],
		filter = config['variant_filter'],
		psrecord = "log/psrecord/joint_genotyping/variant_filtration/variant_filtration.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/variant_filtration/variant_filtration.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/variant_filtration/variant_filtration.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -T VariantFiltration -R {params.ref_genome} -L {params.chr} {params.filter} -V {input.vcf} -o {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule remove_failed:
	input:
		vcf = "data/derived_data/joint_genotyping/variant_filtration/variant_filtration.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/variant_filtration/variant_filtration.{chr}.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['select_variants_gc'],
		nt = config['select_variants_nt'],
		chr = "{chr}",
		java_tmp = "temp/joint_genotyping/remove_failed/{chr}",
		gatk_path = config['gatk_path'],
		psrecord = "log/psrecord/joint_genotyping/remove_failed/remove_failed.{chr}.log"
	output:
		vcf = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -L {params.chr} -T SelectVariants -nt {params.nt} -R {params.ref_genome} -ef --removeUnusedAlternates -env -V {input.vcf} -o {output.vcf}" --log {params.psrecord} --include-children --interval 5
		"""
