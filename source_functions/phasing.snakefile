# snakemake -s source_functions/phasing.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" --until phasing_qc phasing_qc_x phasing_qc_y -p &> log/snakemake_log/joint_genotyping/200710.phasing.log

import os

include: "plink_qc.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

os.makedirs("log/slurm_out/phasing", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/phasing/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['phasing_rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
		expand("data/derived_data/joint_genotyping/phasing/phasing.{autosome}.haps.gz", autosome = list(range(1,30))), "data/derived_data/joint_genotyping/phasing/phasing.X.haps.gz", "data/derived_data/joint_genotyping/phasing/phasing.Y.haps.gz"

rule shapeit_sex_check:
	input:
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.fam"
	params:
		shapeit_module = config['shapeit_module'],
		psrecord = "log/psrecord/joint_genotyping/shapeit_sex_check/shapeit_sex_check.log"
	output:
		log = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.log",
		snp_hh = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.snp.hh",
		ind_hh = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.ind.hh"
	shell:
		"""
		module load {params.shapeit_module}
		psrecord "shapeit check --input-bed {input.bed} {input.bim} {input.fam} --chrX --output-log {output.log}" --log {params.psrecord} --include-children --interval 5
		"""

rule phase_autosomes:
	input:
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.{autosome}.fam",
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
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/plink_qc.X.fam",
		genetic_map = config['genetic_map'],
		log = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.log",
		snp_hh = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.snp.hh",
		ind_hh = "data/derived_data/joint_genotyping/shapeit_sex_check/shapeit_sex_check.ind.hh"
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
		bed = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/plink_qc.Y.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/phasing_qc.Y.fam",
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
