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
		"data/derived_data/joint_genotyping/bovine_demo.guess_ploidy.txt", expand("data/derived_data/joint_genotyping/snp_positions/snp_positions.{chr}.txt", chr = config['chr'])

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
		psrecord = "log/psrecord/joint_genotyping/guess_ploidy/guess_ploidy.log"
	output:
		sex_list = "data/derived_data/joint_genotyping/bovine_demo.guess_ploidy.txt"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools +guess-ploidy -r X:1-133300517 --tag GT {input.x_vcf} > {output.sex_list}" --log {params.psrecord} --include-children --interval 5
		"""
