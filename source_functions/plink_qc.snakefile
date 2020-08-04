# snakemake -s source_functions/plink_qc.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200804.plink_qc.log

import os

# include path is relative to the path of this file
include: "post_process.snakefile"

include: "sex_imputation.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/plink_qc", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/plink_qc/{rules}", rules = config['plink_qc_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['plink_qc_rules']):
    os.makedirs(x, exist_ok = True)

rule plink_qc_all:
	input:
		expand("data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.bed", autosome = list(range(1,30))), "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.bed", "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bed"

rule qc_autosomes:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{autosome}.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		prefix = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter'],
		autosome = "{autosome}"
	output:
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr {params.autosome} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

rule qc_x:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.X.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		pab = config['pab'],
		prefix = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr X -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --from-bp 1 --to-bp {params.pab} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

# set heterozygous calls to missing on the Y
rule qc_y:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.Y.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		prefix = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y",
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr Y --set-hh-missing --geno {params.geno_filter} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --out {params.prefix}
		"""
