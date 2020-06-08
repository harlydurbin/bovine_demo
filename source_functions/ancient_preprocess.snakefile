#snakemake -s source_functions/ancient_preprocess.snakefile --rerun-incomplete --jobs 434 --latency-wait 60 --config --cluster-config source_functions/cluster/ancient_preprocess.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/ancient_preprocess/200522.ancient_preprocess.log

 # 8)  GATK RealignerTargetCreator
 # 9)  GATK IndelRealigner
 # 10) GATK BQSR
 # 11) GATK HaplotypeCaller gVCF mode

#/storage/hpc/group/UMAG/WORKING/hjdzpd/bovine_diversity/data/raw_data/aurochs

# scp -vvv *.bai hjdzpd@lewis.rnet.missouri.edu:/storage/hpc/group/UMAG/WORKING/hjdzpd/bovine_diversity/data/raw_data/aurochs

# scp acad-colab1@acad2.rc.adelaide.edu.au:/localscratch/BisonProjects/Genomes/ARS1_ALN/aurochs/paleomix/*.bai data/raw_data/aurochs

import os

configfile : "source_functions/config/ancient_preprocess.config.yaml"

# Make log directories if they don't exist
for x in expand("log/slurm_out/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/ancient_preprocess", exist_ok = True)

# Make temp directories if they don't exist
for x in expand("temp/ancient_preprocess/{pop}", pop = config['pop']):
    os.makedirs(x, exist_ok = True)

rule all:
	input:
		expand("data/derived_data/ancient_preprocess/{pop}/haplotype_caller.{pop}.{chr}.g.vcf.gz",
		pop = config["pop"],
		chr = config["chr"])

bam_path = {"siberian": "data/raw_data/aurochs/A2494_Bprimigenius_RussiaYeniseiRiver_14082.ARS1_UCD1.2.realigned.bam",
"derbyshire": "data/raw_data/aurochs/CPC98_Bprimigenius_EnglandDerbyshire_5936.ARS1_UCD1.2.realigned.bam"}

# rule to trim reads
# "trimBam will modify the sequences to 'N', and the quality string to '!' unless the optional parameter --clip/-c is specified. If --clip/-c is specified, the ends will be soft clipped instead of modified."

rule trimBam:
	input:
		bam = lambda wildcards: expand("{path}", path = bam_path[wildcards.pop])
	params:
		psrecord = "log/psrecord/ancient_preprocess/trimBam.{pop}.log",
		tmpdir = "temp/ancient_preprocess/{pop}",
		bamUtil_path = config["bamUtil_path"],
		nbases = 5
	output:
		trimmed = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam"
	shell:
	#-c is soft trimming
	# Don't soft trim for now
		"""
		psrecord "{params.bamUtil_path} trimBam {input.bam} {output.trimmed} {params.nbases}" --log {params.psrecord} --include-children --interval 5
		"""

rule index_bam:
	input:
		bam = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam"
	params:
		psrecord = "log/psrecord/ancient_preprocess/index_bam.{pop}.log",
		samtools_module = config['samtools_module']
	output:
		bai = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam.bai"
	shell:
		"""
		module load {params.samtools_module}
		psrecord "samtools index {input.bam}" --log {params.psrecord} --include-children --interval 5
		"""

rule target_creator:
	input:
		indexed_bam = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam",
		bai = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam.bai"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['target_creator_gc'],
		xmx = config['target_creator_xmx'],
		java_tmp = "temp/ancient_preprocess/{pop}",
		gatk_path = config['gatk_path'],
		nt = config['target_creator_nt'],
		chr = "{chr}",
		psrecord = "log/psrecord/ancient_preprocess/target_creator.{pop}.{chr}.log"
	output:
		intervals = "data/derived_data/ancient_preprocess/{pop}/target_creator.{pop}.{chr}.intervals"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -nt {params.nt} -T RealignerTargetCreator -R {params.ref_genome} -L {params.chr} -I {input.indexed_bam} -o {output.intervals}" --log {params.psrecord} --include-children --interval 5
		"""

rule indel_realigner:
	input:
		indexed_bam = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam",
		bai = "data/derived_data/ancient_preprocess/{pop}/trimBam.{pop}.bam.bai",
		intervals = "data/derived_data/ancient_preprocess/{pop}/target_creator.{pop}.{chr}.intervals"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['indel_realigner_gc'],
		xmx = config['indel_realigner_xmx'],
		java_tmp = "temp/ancient_preprocess/{pop}",
		gatk_path = config['gatk_path'],
		chr = "{chr}",
		psrecord = "log/psrecord/ancient_preprocess/indel_realigner.{pop}.{chr}.log"
	output:
		realigned_bam = "data/derived_data/ancient_preprocess/{pop}/indel_realigner.{pop}.{chr}.bam"
	shell:
		"""
		module load {params.java_module}
		module load {params.gatk_module}
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -XX:ParallelGCThreads={params.gc_threads} -Xmx{params.xmx}g -jar {params.gatk_path} -R {params.ref_genome} -I {input.indexed_bam} -T IndelRealigner -targetIntervals {input.intervals} -L {params.chr} -o {output.realigned_bam} --bam_compression 0" --log {params.psrecord} --include-children --interval 5
		"""
rule merge_realigned:
	input:
		realigned_bam = expand("data/derived_data/ancient_preprocess/{{pop}}/indel_realigner.{{pop}}.{chr}.bam", chr = config['chr'])
	params:
		realigned_wildcard = "data/derived_data/ancient_preprocess/{pop}/indel_realigner.{pop}.*.bam",
		nt = config['merge_realigned_nt'],
		psrecord = "log/psrecord/ancient_preprocess/merge_realigned.{pop}.log",
		samtools_module = config['samtools_module']
	output:
		merge_list = "data/derived_data/ancient_preprocess/{pop}/merge_realigned.files.{pop}.list",
		merged_bam = "data/derived_data/ancient_preprocess/{pop}/merge_realigned.{pop}.bam"
	shell:
		"""
		module load {params.samtools_module}
		psrecord "ls {params.realigned_wildcard} | tr '\t' '\n' > {output.merge_list};
		samtools merge -@ {params.nt} -u -f -c -b {output.merge_list} {output.merged_bam};
		samtools index {output.merged_bam}" --log {params.psrecord} --include-children --interval 5
		"""

# Call haplotypes chromosome by chromosome
rule haplotype_caller:
	input:
		merged_bam = "data/derived_data/ancient_preprocess/{pop}/merge_realigned.{pop}.bam"
	params:
		java_module = config['java_module'],
		gatk_module = config['gatk_module'],
		ref_genome = config['ref_genome'],
		gc_threads = config['haplotype_caller_gc'],
		xmx = config['haplotype_caller_xmx'],
		nt = config['haplotype_caller_nt'],
		java_tmp = "temp/ancient_preprocess/{pop}/",
		gatk_path = config['gatk_path'],
		chr = "{chr}",
		psrecord = "log/psrecord/ancient_preprocess/haplotype_caller/haplotype_caller.{pop}.{chr}.log"
	output:
		gvcf = "data/derived_data/ancient_preprocess/{pop}/haplotype_caller.{pop}.{chr}.g.vcf.gz",
		tbi = "data/derived_data/ancient_preprocess/{pop}/haplotype_caller.{pop}.{chr}.g.vcf.gz.tbi"
	shell:
		"""
		psrecord "java -Djava.io.tmpdir={params.java_tmp} -Xmx{params.xmx}g -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_path} -nct {params.nt} -ERC GVCF -T HaplotypeCaller -R {params.ref_genome} -L {params.chr} -I {input.merged_bam} --pcr_indel_model NONE --heterozygosity 0.0033 -o {output.gvcf}" --log {params.psrecord} --include-children --interval 5
		"""
