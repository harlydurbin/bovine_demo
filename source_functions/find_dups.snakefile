# snakemake -s source_functions/find_dups.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/find_dups.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200629.find_dups.log

import os

configfile: "source_functions/config/find_dups.config.yaml"

os.makedirs("log/slurm_out/find_dups", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/find_dups/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("temp/find_dups", exist_ok = True)

rule all:
	input:
	 	"data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.con"

# rule index_map:
# 	input:
# 		map = config['map']
# 	output:
# 		regions_file = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.list"
# 	shell:
# 		"""
# 		awk '{{print $1":"$2}}' {input.map} > {output.regions_file}
		# """

rule index_map:
	input:
		map = config['map']
	params:
		bcftools_module = config['bcftools_module']
	output:
		regions_file = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.tsv.gz",
		regions_file_tbi = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.tsv.gz.tbi"
	shell:
		"""
		module load {params.bcftools_module}
		bgzip -c {input.map} > {output.regions_file}
		tabix -s1 -b2 -e2 {output.regions_file}
		"""

rule extract_850K:
	input:
		regions_file = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.tsv.gz",
		regions_file_tbi = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.tsv.gz.tbi",
		full_file = config['full_file'],
	params:
		psrecord = "log/psrecord/joint_genotyping/extract_850K/extract_850K.log",
		bcftools_module = config['bcftools_module']
	output:
		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools view --regions-file {input.regions_file} -O z -o {output.subset_file} {input.full_file}" --log {params.psrecord} --include-children --interval 5
		"""

# rule extract_850K:
# 	input:
# 		full_file = config['full_file'],
# 		regions_file = "data/derived_data/joint_genotyping/find_dups/regions_file.850K.list"
# 	params:
# 		psrecord = "log/psrecord/joint_genotyping/extract_850K/extract_850K.log",
# 		java_module = config['java_module'],
# 		gatk_path = config['gatk_path'],
# 		java_tmp = "temp/find_dups",
# 		gc_threads = config['extract_850K_gc'],
# 		nt = config['extract_850K_nt'],
# 		ref_genome = config['ref_genome']
# 	output:
# 		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz",
# 		tbi = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz.tbi"
# 	shell:
# 		"""
# 		module load {params.java_module}
# 		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -L {input.regions_file} -T SelectVariants -nt {params.nt} -R {params.ref_genome} -V {input.full_file} -o {output.subset_file}" --log {params.psrecord} --include-children --interval 5
# 		"""


rule index_850K:
	input:
		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz"
	params:
		bcftools_module = config['bcftools_module'],
		nt = config['index_850K_nt'],
		psrecord = "log/psrecord/joint_genotyping/index_850K/index_850K.log"
	output:
		tbi = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz.tbi"
	shell:
		"""
		module load {params.bcftools_module}
		psrecord "bcftools index --tbi -f --threads {params.nt} {input.subset_file}" --log {params.psrecord} --include-children --interval 5
		"""

rule dups_plink:
	input:
		subset_file = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz",
		tbi = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.vcf.gz.tbi"
	params:
		prefix = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K",
		nt = config['plink_nt'],
		psrecord = "log/psrecord/joint_genotyping/dups_plink/dups_plink.log"
	output:
		bed = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.bed"
	shell:
		"""
		module load plink
		psrecord "plink --vcf {input.subset_file} --make-bed --double-id --cow --threads {params.nt} --out {params.prefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule find_dups:
	input:
		bed = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.bed"
	params:
		king_path = config['king_path'],
		prefix = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K",
		nt = config['king_nt'],
		psrecord = "log/psrecord/joint_genotyping/find_dups/find_dups.log"
	output:
		dups = "data/derived_data/joint_genotyping/find_dups/bovine_demo.850K.con"
	shell:
		"""
		psrecord "{params.king_path} -b {input.bed} --duplicate --prefix {params.prefix} --sexchr 30 --cpus {params.nt}" --log {params.psrecord} --include-children --interval 5
		"""
