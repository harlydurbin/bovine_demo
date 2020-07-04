# snakemake -s source_functions/find_dups.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/genotyping_qc_phasing.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200704.find_dups.log

# include path is relative to the path of this file
include: "joint_genotyping.snakefile"

configfile: "source_functions/config/genotyping_qc_phasing.config.yaml"

import os

os.makedirs("log/slurm_out/find_dups", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/find_dups/{rules}", rules = config['find_dups_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['find_dups_rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("temp/find_dups", exist_ok = True)

rule find_dups_all:
	input:
	 	"data/derived_data/joint_genotyping/find_dups/find_dups.con"

rule targets_file:
	input:
		map = config['850K_map']
	params:
		chr = "{chr}"
	output:
		targets_file = "data/derived_data/joint_genotyping/find_dups/{chr}.targets"
	shell:
		"""
		grep "^{params.chr}" {input.map} | awk '{{print $1"\t"$2}}' > {output.targets_file}
		"""

rule extract_850K:
	input:
		targets_file = "data/derived_data/joint_genotyping/find_dups/{chr}.targets",
		vcf = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/remove_failed/remove_failed.{chr}.vcf.gz.tbi"
	params:
		psrecord = "log/psrecord/joint_genotyping/extract_850K/extract_850K.{chr}.log",
		bcftools_module = config['bcftools_module'],
		chr = "{chr}"
	output:
		subset_file = "data/derived_data/joint_genotyping/find_dups/extract_850K.{chr}.vcf.gz"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools view --regions {params.chr} --targets-file {input.targets_file} -O z -o {output.subset_file} {input.vcf}" --log {params.psrecord} --include-children --interval 5
		"""

rule index_850K:
	input:
		subset_file = "data/derived_data/joint_genotyping/find_dups/extract_850K.{chr}.vcf.gz"
	params:
		bcftools_module = config['bcftools_module'],
		nt = config['index_850K_nt'],
		psrecord = "log/psrecord/joint_genotyping/index_850K/index_850K.{chr}.log"
	output:
		tbi = "data/derived_data/joint_genotyping/find_dups/extract_850K.{chr}.vcf.gz.tbi"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools index --tbi -f --threads {params.nt} {input.subset_file}" --log {params.psrecord} --include-children --interval 5
		"""

rule make_bed:
	input:
		subset_file = "data/derived_data/joint_genotyping/find_dups/extract_850K.{chr}.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/find_dups/extract_850K.{chr}.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/find_dups/make_bed.{chr}",
		nt = config['make_bed_nt'],
		psrecord = "log/psrecord/joint_genotyping/make_bed/make_bed.{chr}.log"
	output:
		chr_bed = "data/derived_data/joint_genotyping/find_dups/make_bed.{chr}.bed"
	# Need -set-missing-var-ids in order to avoid SNP id errors!
	shell:
		"""
		module load plink
		psrecord "plink --vcf {input.subset_file} --make-bed --double-id --cow -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --out {params.prefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule merge_list:
	input:
	# Don't worry about sex chromosomes
		chr_bed = expand("data/derived_data/joint_genotyping/find_dups/make_bed.{chr}.bed", chr = list(range(1,30)))
	params:
		prefixes = lambda wildcards: expand("data/derived_data/joint_genotyping/find_dups/make_bed.{chr}\n", chr = list(range(1,30)))
	output:
		merge_list = "data/derived_data/joint_genotyping/find_dups/850K_merge_list.txt"
	shell:
		"echo -e '{params.prefixes}' | sed 's/^ *//g' > {output.merge_list}"

rule merge_bed:
	input:
		chr_bed = expand("data/derived_data/joint_genotyping/find_dups/make_bed.{chr}.bed", chr =  list(range(1,30))),
		merge_list = "data/derived_data/joint_genotyping/find_dups/850K_merge_list.txt"
	params:
		prefix = "data/derived_data/joint_genotyping/find_dups/merge_bed",
		nt = config['merge_bed_nt'],
		psrecord = "log/psrecord/joint_genotyping/merge_bed/merge_bed.log"
	output:
		bed = "data/derived_data/joint_genotyping/find_dups/merge_bed.bed"
	shell:
		"""
		module load plink
		psrecord "plink --merge-list {input.merge_list} --make-bed --double-id --cow --threads {params.nt} --out {params.prefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule find_dups:
	input:
		bed = "data/derived_data/joint_genotyping/find_dups/merge_bed.bed"
	params:
		king_path = config['king_path'],
		prefix = "data/derived_data/joint_genotyping/find_dups/find_dups",
		nt = config['king_nt'],
		psrecord = "log/psrecord/joint_genotyping/find_dups/find_dups.log"
	output:
		dups = "data/derived_data/joint_genotyping/find_dups/find_dups.con"
	shell:
		"""
		psrecord "{params.king_path} -b {input.bed} --duplicate --prefix {params.prefix} --sexchr 30 --cpus {params.nt}" --log {params.psrecord} --include-children --interval 5
		"""
