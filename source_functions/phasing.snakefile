# snakemake -s source_functions/phasing.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200708.phasing.log

import os

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/phasing", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/phasing/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("temp/phasing", exist_ok = True)

rule all:
	input:
		expand("data/derived_data/joint_genotyping/snp_positions/snp_positions.{chr}.txt", chr = config['chr']), expand("data/derived_data/joint_genotyping/impute_sex/bovine_demo.guess_ploidy.{tag}.txt", tag = config['tag']), expand("data/derived_data/joint_genotyping/phasing/phasing.{autosome}.vcf.gz", autosome = list(range(1,30)) "Y")), "data/derived_data/joint_genotyping/phasing/phasing.X.vcf.gz"

rule snp_positions:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{chr}.vcf.gz"
	params:
		bcftools_module = config['bcftools_module']
	output:
		pos_list = "data/derived_data/joint_genotyping/snp_positions/snp_positions.{chr}.txt"
	shell:
		"""
		module load {params.bcftools_module}
		bcftools query -f "%CHROM\\t%POS\\n" {input.vcf} >  {output.pos_list}
		"""

rule guess_ploidy:
	input:
		x_vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.X.vcf.gz"
	params:
		bcftools_module = config['bcftools_module'],
		psrecord = "log/psrecord/joint_genotyping/guess_ploidy/guess_ploidy.{tag}.log",
		tag = "{tag}"
	output:
		sex_list = "data/derived_data/joint_genotyping/impute_sex/bovine_demo.guess_ploidy.{tag}.txt"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools +guess-ploidy -r X:1-133300517 --tag {params.tag} {input.x_vcf} > {output.sex_list}" --log {params.psrecord} --include-children --interval 5
		"""

rule plink_impute_sex:
	input:
		x_vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.X.vcf.gz"
	params:
		plink_module = config['plink_module'],
	output:
		fam = "data/derived_data/joint_genotyping/impute_sex/impute_sex.fam",
		sexcheck = "data/derived_data/joint_genotyping/impute_sex/impute_sex.sexcheck"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr X -set-missing-var-ids @:#\$1\$2 --threads 12 --from-bp 1 --to-bp 133300518 --out data/derived_data/joint_genotyping/impute_sex/remov_par
		plink --bfile data/derived_data/joint_genotyping/impute_sex/remov_par --double-id --cow --threads 12 --impute-sex --make-bed --out data/derived_data/joint_genotyping/impute_sex/impute_sex
		"""

rule phasing_autosomes:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{autosome}.vcf.gz",
		genetic_map = config['genetic_map']
	output:
		phased_vcf = "data/derived_data/joint_genotyping/phasing/phasing.{autosome}.vcf.gz"

rule phasing_x:
	input:
		x_vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.X.vcf.gz",
		genetic_map = config['genetic_map'],
		sex_list = config['sex_list']
	output:
		phased_vcf = "data/derived_data/joint_genotyping/phasing/phasing.X.vcf.gz"
