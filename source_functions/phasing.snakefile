# snakemake -s source_functions/phasing.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" --until phasing_qc phasing_qc_x phasing_qc_y -p &> log/snakemake_log/joint_genotyping/200710.phasing.log

import os

# include path is relative to the path of this file
include: "post_process.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/phasing", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/phasing/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
		expand("data/derived_data/joint_genotyping/snp_positions/snp_positions.{chr}.txt", chr = config['chr']), expand("data/derived_data/joint_genotyping/impute_sex/bovine_demo.guess_ploidy.{tag}.txt", tag = config['tag']), expand("data/derived_data/joint_genotyping/phasing/phasing.{autosome}.haps.gz", autosome = list(range(1,30))), "data/derived_data/joint_genotyping/phasing/phasing.X.haps.gz", "data/derived_data/joint_genotyping/phasing/phasing.Y.haps.gz"

#### PHASING ####

rule phasing_qc:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.{autosome}.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		prefix = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter'],
		autosome = "{autosome}"
	output:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr {params.autosome} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

rule phasing_qc_x:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.X.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		pab = config['pab'],
		prefix = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr X -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --from-bp 1 --to-bp {params.pab} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

# set heterozygous calls to missing on the Y
rule phasing_qc_y:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/remove_samples.Y.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		prefix = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y",
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr Y --set-hh-missing --geno {params.geno_filter} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --out {params.prefix}
		"""

rule shapeit_sex_check:
	input:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.fam"
	params:
		shapeit_module = config['shapeit_module'],
		psrecord = "log/psrecord/joint_genotyping/shapeit_sex_check/shapeit_sex_check.log"
	output:
		log = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.log",
		snp_hh = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.snp.hh",
		ind_hh = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.ind.hh"
	shell:
		"""
		module load {params.shapeit_module}
		psrecord "shapeit check --input-bed {input.bed} {input.bim} {input.fam} --chrX --output-log {output.log}" --log {params.psrecord} --include-children --interval 5
		"""

rule phase_autosomes:
	input:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{autosome}.fam",
		genetic_map = config['genetic_map']
	params:
		shapeit_module = config['shapeit_module'],
		psrecord = "log/psrecord/joint_genotyping/phase_autosomes/phase_autosomes.{autosome}.log",
		nt = config['shapeit_nt']
	output:
		haps = "data/derived_data/joint_genotyping/phasing/phasing.{autosome}.haps.gz",
		sample = "data/derived_data/joint_genotyping/phasing/phasing.{autosome}.sample"
	shell:
		"""
		module load {params.shapeit_module}
		psrecord "shapeit --input-bed {input.bed} {input.bim} {input.fam} -M {input.genetic_map} -O {output.haps} {output.sample} --thread {params.nt}" --log {params.psrecord} --include-children --interval 5
		"""

rule phase_x:
	input:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.X.fam",
		genetic_map = config['genetic_map'],
		log = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.log",
		snp_hh = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.snp.hh",
		ind_hh = "data/derived_data/joint_genotyping/phasing_qc/shapeit_sex_check.ind.hh"
	params:
		shapeit_module = config['shapeit_module'],
		psrecord = "log/psrecord/joint_genotyping/phase_x/phase_x.log",
		nt = config['shapeit_nt']
	output:
		haps = "data/derived_data/joint_genotyping/phasing/phasing.X.haps.gz",
		sample = "data/derived_data/joint_genotyping/phasing/phasing.X.sample"
	shell:
		"""
		module load {params.shapeit_module}
		psrecord "shapeit --input-bed {input.bed} {input.bim} {input.fam} -M {input.genetic_map} -O {output.haps} {output.sample} --chrX --thread {params.nt}" --log {params.psrecord} --include-children --interval 5
		"""

rule phase_y:
	input:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.Y.fam",
		genetic_map = config['genetic_map']
	params:
		shapeit_module = config['shapeit_module'],
		psrecord = "log/psrecord/joint_genotyping/phase_y/phase_y.log",
		nt = config['shapeit_nt']
	output:
		haps = "data/derived_data/joint_genotyping/phasing/phasing.Y.haps.gz",
		sample = "data/derived_data/joint_genotyping/phasing/phasing.Y.sample"
	shell:
		"""
		module load {params.shapeit_module}
		psrecord "shapeit --input-bed {input.bed} {input.bim} {input.fam} -M {input.genetic_map} -O {output.haps} {output.sample} --thread {params.nt}" --log {params.psrecord} --include-children --interval 5
		"""
