# snakemake -s source_functions/sex_imputation.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200728.sex_imputation.log

import os

# include path is relative to the path of this file
# include: "post_process.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/sex_imputation", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/sex_imputation/{rules}", rules = config['sex_imputation_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['sex_imputation_rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
		expand("data/derived_data/joint_genotyping/sex_imputation/bovine_demo.guess_ploidy.{tag}.txt", tag = config['tag']), "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.impute_sex.sexcheck", "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.impute_sex.fam"

rule guess_ploidy:
	input:
		x_vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.X.vcf.gz"
	params:
		bcftools_module = config['bcftools_module'],
		psrecord = "log/psrecord/joint_genotyping/guess_ploidy/guess_ploidy.{tag}.log",
		tag = "{tag}",
		pab = config['pab']
	output:
		sex_list = "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.guess_ploidy.{tag}.txt"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools +guess-ploidy -r X:1-{params.pab} --tag {params.tag} {input.x_vcf} > {output.sex_list}" --log {params.psrecord} --include-children --interval 5
		"""

rule plink_impute_sex:
	input:
		x_vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.X.vcf.gz"
	params:
		plink_module = config['plink_module'],
		pab = config['pab'],
		remove_par_prefix = "data/derived_data/joint_genotyping/sex_imputation/remove_par",
		impute_sex_prefix = "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.impute_sex"
	output:
		fam = "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.impute_sex.fam",
		sexcheck = "data/derived_data/joint_genotyping/sex_imputation/bovine_demo.impute_sex.sexcheck"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.x_vcf} --make-bed --double-id --cow --chr X -set-missing-var-ids @:#\$1\$2 --threads 12 --from-bp 1 --to-bp {params.pab} --out {params.remove_par_prefix}
		plink --bfile {params.remove_par_prefix} --double-id --cow --threads 12 --impute-sex --make-bed --out {params.impute_sex_prefix}
		"""
