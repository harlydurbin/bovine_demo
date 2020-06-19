# snakemake -s source_functions/find_dups.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/find_dups.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200619.find_dups.log

import os

configfile: "source_functions/config/find_dups.config.yaml"

# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
	 	"data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.con"

rule index_map:
	input:
		map = config['map']
	params:
		bcftools_module = config['bcftools_module']
	output:
		map_gz = config['map'] + ".gz",
		map_index = config['map'] + ".gz.tbi"
	shell:
		"""
		module load {params.bcftools_module}
		bgzip -c {input.map} > {output.map_gz}
		tabix -s1 -b2 -e2 {output.map_gz}
		"""

rule extract_variants:
	input:
		full_file = config['full_file'],
		map_gz = config['map'] + ".gz",
		map_index = config['map'] + ".gz.tbi"
	params:
		bcftools_module = config['bcftools_module']
	output:
		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz"
	shell:
		"""
		module load {params.bcftools_module}
		bcftools view --regions-file {input.map_gz} -O z -o {output.subset_file} {input.full_file}
		tabix {output.subset_file}
		"""

rule dups_plink:
	input:
		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz"
	params:
		prefix = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K",
		nt = config['plink_nt']
	output:
		bed = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.bed"
	shell:
		"""
		module load plink
		plink --bcf {input.subset_file} --make-bed --double-id --cow --threads {params.nt} --out {params.prefix}
		"""

rule find_dups:
	input:
		bed = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.bed"
	params:
		king_path = config['king_path'],
		prefix = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K",
		nt = config['king_nt']
	output:
		dups = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.con"
	shell:
		"{params.king_path} -b {input.bed} --duplicate --prefix {params.prefix} --sexchr 30 --cpus {params.nt}"
