import os

configfile: "source_functions/config/filter_eval.config.yaml"

# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("temp/filter_eval/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
	 	"data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.table", "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.ldepth.mean", "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.con", "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.gf.table"

# Create a chr28 vcf with snps & indels to evaluate how many indels are within
# 5bp of an inel, filtering values
rule biallelic_28:
	input:
		vcf = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/genotype_gvcfs/genotype_gvcfs.28.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		nt = config['biallelic_28_nt'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/filter_eval/biallelic_28"
	output:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz.tbi"
	shell:
		"""
		module load {params.java_module}
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 {params.gatk_path} -nt {params.nt} -T SelectVariants -R {params.ref_genome} -L 28 -V {input.vcf} --restrictAllelesTo BIALLELIC -o {output.vcf}
		"""

rule plink_28:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28",
		nt = config['plink_nt']
	output:
		bed = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.bed"
	shell:
		"""
		module load plink
		plink --vcf {input.vcf} --make-bed --double-id --cow --threads {params.nt} --out {params.prefix}
		"""

# vcftools --gzvcf data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz --get-INFO GQ --out data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28

rule depth_28:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28"
	output:
		depth = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.ldepth.mean"
	shell:
		"""
		module load vcftools
		vcftools --gzvcf {input.vcf} --site-mean-depth --out {params.prefix}
		"""

rule table_28:
	input:
		vcf = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.vcf.gz.tbi"
	params:
		java_module = config['java_module'],
		ref_genome = config['ref_genome'],
		gatk_path = config['gatk_path'],
		java_tmp = "temp/joint_genotyping/table_28"
	output:
		table = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.table"
	shell:
		"""
		module load {params.java_module}
		java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads=2 -jar {params.gatk_path} -R {params.ref_genome} -L 28 -T VariantsToTable -V {input.vcf} -F POS -F TYPE -F TRANSITION -F QD -F FS -F MQ -F ReadPosRankSum -F MQRankSum -F NO-CALL -F N-CALLED -F VAR -o {output.table}
		"""

rule find_dups:
	input:
		bed = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.bed"
	params:
		king_path = config['king_path'],
		prefix = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28",
		nt = config['king_nt']
	output:
		dups = "data/derived_data/joint_genotyping/filter_eval/filter_eval.with_indels.28.con"
	shell:
		"{params.king_path} -b {input.bed} --duplicate --prefix {params.prefix} --sexchr 30 --cpus {params.nt}"
