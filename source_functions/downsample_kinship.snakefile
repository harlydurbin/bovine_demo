#snakemake -s source_functions/downsample_kinship.snakefile -j 250 --rerun-incomplete --latency-wait 60 --config --cluster-config source_functions/cluster/downsample_kinship.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -n {cluster.n} --mem {cluster.mem} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/190314.downsample_kinship.log

#snakemake -s source_functions/downsample_kinship.snakefile --forceall --rulegraph | dot -Tpdf > data/derived_data/sample_selection/downsample_kinship.rulegraph.pdf

#smartpca -p {input.parfile} > {output.log}

#smartpca -p data/derived_data/sample_selection/ds_smartpca/angus.par data/derived_data/sample_selection/ds_smartpca/angus

#downsampling based on kinship
configfile : "source_functions/config/downsample_kinship.config.json"

rule target:
	input:
		targ = expand("data/derived_data/sample_selection/ds_grm/{pop}.ds_grm.sXX.txt",
		pop = config['pops']),
		targ2 = "data/derived_data/sample_selection/ds_grm_big/ds_grm_big.sXX.txt"

#dictionary of how many cohorts per population
#things screw up when these are just numbers instead of prefixed by "cohort"
cohort_dict = { "angus": ["cohort1", "cohort2", "cohort3", "cohort4", "cohort5"],
"brahman":["cohort1"],
"charolais": ["cohort1", "cohort2"],
"hereford": ["cohort1", "cohort2"],
"holstein":["cohort1", "cohort2", "cohort3", "cohort4", "cohort5"],
"jersey": ["cohort1", "cohort2"],
"limousin": ["cohort1"],
"simmental": ["cohort1", "cohort2", "cohort3"]
}

big_list = ["angus.cohort1",
"angus.cohort2",
"angus.cohort3",
"angus.cohort4",
"angus.cohort5",
"brahman.cohort1",
"charolais.cohort1",
"charolais.cohort2",
"hereford.cohort1",
"hereford.cohort2",
"holstein.cohort1",
"holstein.cohort2",
"holstein.cohort3",
"holstein.cohort4",
"holstein.cohort5",
"jersey.cohort1",
"jersey.cohort2",
"limousin.cohort1",
"simmental.cohort1",
"simmental.cohort2",
"simmental.cohort3"]


rule ds_combine_gvcf:
	input:
		gvcflist = "data/derived_data/sample_selection/downsample_cohorts/{pop}.{cohort}.list"
	benchmark:
		"benchmarks/{rule}/{pop}.{cohort}.{rule}.benchmark"
	# threads:
	# 	config['combine_gvcf_nt'] # can't parallelize combinegvcf
	params:
		tmpdir = "temp/{rule}/{pop}.{cohort}/",
		cohort = "{cohort}",
		ref_genome = config['ref_genome'],
		gatk_version = config['gatk_version'],
		gc_threads = config["combine_gvcf_gc"],
		# nt = config["combine_gvcf_nt"], # can't parallelize combinegvcf
		# xmx = config["combine_gvcf_xmx"],
		gatklog = "log/gatklog/{rule}/{pop}.{cohort}.{rule}.log"
	log:
		"log/snakemake_log/{rule}/{pop}.{cohort}.{rule}.snakelog",
	output:
		gvcf = "data/derived_data/sample_selection/{rule}/{pop}.{cohort}.{rule}.g.vcf.gz",
		tbi = "data/derived_data/sample_selection/{rule}/{pop}.{cohort}.{rule}.g.vcf.gz.tbi"
	shell:
		"java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T CombineGVCFs -R  {params.ref_genome} -V {input.gvcflist} --log_to_file {params.gatklog} -o {output.gvcf}"

rule ds_genotype_gvcf:
	input:
		combine_gvcf = lambda wildcards: expand('data/derived_data/sample_selection/ds_combine_gvcf/{pop}.{cohort}.ds_combine_gvcf.g.vcf.gz', pop = wildcards.pop, cohort = cohort_dict[wildcards.pop]),
		combine_tbi = lambda wildcards: expand('data/derived_data/sample_selection/ds_combine_gvcf/{pop}.{cohort}.ds_combine_gvcf.g.vcf.gz.tbi', pop = wildcards.pop, cohort = cohort_dict[wildcards.pop])
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/{pop}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["genotype_gvcf_gc"],
		nt = config['genotype_gvcf_nt'],
		# xmx = config["genotype_gvcf_xmx"],
		combine_gvcf = lambda wildcards: expand("-V data/derived_data/sample_selection/ds_combine_gvcf/{pop}.{cohort}.ds_combine_gvcf.g.vcf.gz", pop = wildcards.pop, cohort = cohort_dict[wildcards.pop]),
		gatklog = "log/gatklog/{rule}/{pop}.{rule}.log"
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		vcf = "data/derived_data/sample_selection/{rule}/{pop}.ds_genotype_gvcf.vcf.gz",
		tbi = "data/derived_data/sample_selection/{rule}/{pop}.ds_genotype_gvcf.vcf.gz.tbi"
	shell:
		"java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -nt {params.nt} -T GenotypeGVCFs -R {params.ref_genome} -L 28 {params.combine_gvcf} --useNewAFCalculator --heterozygosity 0.0033 --standard_min_confidence_threshold_for_calling 10 --log_to_file {params.gatklog} -o {output.vcf}"

#GenotypeGVCFs all together so I can make a lorg GRM
rule ds_genotype_gvcf_big:
	input:
		combine_gvcf = lambda wildcards: expand('data/derived_data/sample_selection/ds_combine_gvcf/{prefix}.ds_combine_gvcf.g.vcf.gz', prefix = big_list),
		combine_gvcf_tbi = lambda wildcards: expand('data/derived_data/sample_selection/ds_combine_gvcf/{prefix}.ds_combine_gvcf.g.vcf.gz.tbi', prefix = big_list)
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	# threads:
	# 	config['combine_gvcf_nt'] # can't parallelize combinegvcf
	params:
		tmpdir = "temp/{rule}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["genotype_gvcf_big_gc"],
		nt = config['genotype_gvcf_big_nt'],
		# xmx = config["genotype_gvcf_xmx"],
		combine_gvcf = lambda wildcards: expand('-V data/derived_data/sample_selection/ds_combine_gvcf/{prefix}.ds_combine_gvcf.g.vcf.gz', prefix = big_list),
		gatklog = "log/gatklog/{rule}/{rule}.log"
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog",
	output:
		vcf = "data/derived_data/sample_selection/{rule}/ds_genotype_gvcf_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/{rule}/ds_genotype_gvcf_big.vcf.gz.tbi"
	shell:
		"java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -nt {params.nt} -T GenotypeGVCFs -R {params.ref_genome} -L 28 {params.combine_gvcf} --useNewAFCalculator --heterozygosity 0.0033 --standard_min_confidence_threshold_for_calling 10 --log_to_file {params.gatklog} -o {output.vcf}"

rule ds_remove_indels:
	input:
		vcf = "data/derived_data/sample_selection/ds_genotype_gvcf/{pop}.ds_genotype_gvcf.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_genotype_gvcf/{pop}.ds_genotype_gvcf.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/{pop}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["remove_indels_gc"],
		gatklog = "log/gatklog/{rule}/{pop}.{rule}.log"
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		vcf = temp("data/derived_data/sample_selection/{rule}/{pop}.ds_remove_indels.vcf.gz"),
		tbi = temp("data/derived_data/sample_selection/{rule}/{pop}.ds_remove_indels.vcf.gz.tbi")
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T SelectVariants -R {params.ref_genome} -selectType SNP -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_remove_indels_big:
	input:
		vcf = "data/derived_data/sample_selection/ds_genotype_gvcf_big/ds_genotype_gvcf_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_genotype_gvcf_big/ds_genotype_gvcf_big.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["remove_indels_big_gc"],
		gatklog = "log/gatklog/{rule}/{rule}.log"
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog"
	output:
		vcf = temp("data/derived_data/sample_selection/{rule}/ds_remove_indels_big.vcf.gz"),
		tbi = temp("data/derived_data/sample_selection/{rule}/ds_remove_indels_big.vcf.gz.tbi")
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T SelectVariants -R {params.ref_genome} -selectType SNP -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_variant_qc:
	input:
		vcf = "data/derived_data/sample_selection/ds_remove_indels/{pop}.ds_remove_indels.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_remove_indels/{pop}.ds_remove_indels.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/{pop}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["variant_qc_gc"],
		gatklog = "log/gatklog/{rule}/{pop}.{rule}.log",
		filter = config["filter"],
		filtername = config["filtername"]
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		vcf = temp("data/derived_data/sample_selection/{rule}/{pop}.ds_variant_qc.vcf.gz"),
		tbi = temp("data/derived_data/sample_selection/{rule}/{pop}.ds_variant_qc.vcf.gz.tbi")
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T VariantFiltration -R {params.ref_genome} --filterExpression {params.filter} --filterName {params.filtername} -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_variant_qc_big:
	input:
		vcf = "data/derived_data/sample_selection/ds_remove_indels_big/ds_remove_indels_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_remove_indels_big/ds_remove_indels_big.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/",
		ref_genome = config["ref_genome"],
		gatk_version = config['gatk_version'],
		gc_threads = config["variant_qc_big_gc"],
		gatklog = "log/gatklog/{rule}/{rule}.log",
		filter = config["filter"],
		filtername = config["filtername"]
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog"
	output:
		vcf = temp("data/derived_data/sample_selection/{rule}/ds_variant_qc_big.vcf.gz"),
		tbi = temp("data/derived_data/sample_selection/{rule}/ds_variant_qc_big.vcf.gz.tbi")
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T VariantFiltration -R {params.ref_genome} --filterExpression {params.filter} --filterName {params.filtername} -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_remove_failed:
	input:
		vcf = "data/derived_data/sample_selection/ds_variant_qc/{pop}.ds_variant_qc.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_variant_qc/{pop}.ds_variant_qc.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/{pop}/",
		ref_genome = config["ref_genome"],
		gatk_version = config["gatk_version"],
		gc_threads = config["remove_failed_gc"],
		gatklog = "log/gatklog/{rule}/{pop}.{rule}.log"
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		vcf = "data/derived_data/sample_selection/{rule}/{pop}.ds_remove_failed.vcf.gz",
		tbi = "data/derived_data/sample_selection/{rule}/{pop}.ds_remove_failed.vcf.gz.tbi"
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T SelectVariants -R {params.ref_genome} -ef -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_remove_failed_big:
	input:
		vcf = "data/derived_data/sample_selection/ds_variant_qc_big/ds_variant_qc_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_variant_qc_big/ds_variant_qc_big.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	params:
		tmpdir = "temp/{rule}/",
		ref_genome = config["ref_genome"],
		gatk_version = config["gatk_version"],
		gc_threads = config["remove_failed_big_gc"],
		gatklog = "log/gatklog/{rule}/{rule}.log"
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog"
	output:
		vcf = "data/derived_data/sample_selection/{rule}/ds_remove_failed_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/{rule}/ds_remove_failed_big.vcf.gz.tbi"
	shell:
		"(java -Djava.io.tmpdir={params.tmpdir} -XX:ParallelGCThreads={params.gc_threads} -jar {params.gatk_version} -T SelectVariants -R {params.ref_genome} -ef -V {input.vcf} -log {params.gatklog} -o {output.vcf}) > {log}"

rule ds_plink:
	input:
		vcf = "data/derived_data/sample_selection/ds_remove_failed/{pop}.ds_remove_failed.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_remove_failed/{pop}.ds_remove_failed.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		path = config["plink_path"],
		oprefix = "data/derived_data/sample_selection/{rule}/{pop}.{rule}",
		nt = config["plink_nt"]
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		bed = "data/derived_data/sample_selection/{rule}/{pop}.{rule}.bed",
		bim = "data/derived_data/sample_selection/{rule}/{pop}.{rule}.bim",
		fam = "data/derived_data/sample_selection/{rule}/{pop}.{rule}.fam"
	shell:
		"{params.path} --vcf {input.vcf} --cow --threads {params.nt} --double-id -set-missing-var-ids @:#\$1\$2 --make-bed --out {params.oprefix}; sed -i -e 's/-9/69/g' {output.fam}"

rule ds_plink_big:
	input:
		vcf = "data/derived_data/sample_selection/ds_remove_failed_big/ds_remove_failed_big.vcf.gz",
		tbi = "data/derived_data/sample_selection/ds_remove_failed_big/ds_remove_failed_big.vcf.gz.tbi"
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	params:
		path = config["plink_path"],
		oprefix = "data/derived_data/sample_selection/{rule}/{rule}",
		nt = config["plink_big_nt"]
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog"
	output:
		bed = "data/derived_data/sample_selection/{rule}/{rule}.bed",
		bim = "data/derived_data/sample_selection/{rule}/{rule}.bim",
		fam = "data/derived_data/sample_selection/{rule}/{rule}.fam"
		#have to replace -9s in plink fam file with something for GEMMA
	shell:
		"{params.path} --vcf {input.vcf} --cow --threads {params.nt} --double-id -set-missing-var-ids @:#\$1\$2 --make-bed --out {params.oprefix}; sed -i -e 's/-9/69/g' {output.fam}"

rule ds_grm:
	input:
		bed = "data/derived_data/sample_selection/ds_plink/{pop}.ds_plink.bed",
		bim = "data/derived_data/sample_selection/ds_plink/{pop}.ds_plink.bim",
		fam = "data/derived_data/sample_selection/ds_plink/{pop}.ds_plink.fam"
	benchmark:
		"benchmarks/{rule}/{pop}.{rule}.benchmark"
	params:
		path = config["gemma_path"],
		inprefix = "data/derived_data/sample_selection/ds_plink/{pop}.ds_plink",
		oprefix = "{pop}.{rule}"
	log:
		"log/snakemake_log/{rule}/{pop}.{rule}.snakelog"
	output:
		log = "data/derived_data/sample_selection/ds_grm/{pop}.{rule}.log.txt",
		sXX = "data/derived_data/sample_selection/ds_grm/{pop}.{rule}.sXX.txt"
	shell:
		"{params.path} -bfile {params.inprefix} -gk 2 -o {params.oprefix}; mv output/{params.oprefix}* data/derived_data/sample_selection/ds_grm/"

rule ds_grm_big:
	input:
		bed = "data/derived_data/sample_selection/ds_plink_big/ds_plink_big.bed",
		bim = "data/derived_data/sample_selection/ds_plink_big/ds_plink_big.bim",
		fam = "data/derived_data/sample_selection/ds_plink_big/ds_plink_big.fam"
	benchmark:
		"benchmarks/{rule}/{rule}.benchmark"
	params:
		path = config["gemma_path"],
		inprefix = "data/derived_data/sample_selection/ds_plink_big/ds_plink_big",
		oprefix = "{rule}"
	log:
		"log/snakemake_log/{rule}/{rule}.snakelog"
	output:
		log = "data/derived_data/sample_selection/{rule}/{rule}.log.txt",
		sXX = "data/derived_data/sample_selection/{rule}/{rule}.sXX.txt"
	shell:
		"{params.path} -bfile {params.inprefix} -gk 2 -o {params.oprefix}; mv output/{params.oprefix}* data/derived_data/sample_selection/ds_grm_big/"

# rule rel_cutoff:
# /cluster/software/plink/plink-1.90b/plink --bfile data/derived_data/sample_selection/ds_plink/{pop}.ds_plink --cow --threads 4 --double-id --rel-cutoff 0.12 --out data/derived_data/sample_selection/rel_cutoff/{pop}.rel_cutoff 

# rule ds_relatedness2:
# 	input:
# 		vcf = "data/derived_data/sample_selection/ds_remove_failed/{pop}.ds_remove_failed.vcf.gz",
# 		tbi = "data/derived_data/sample_selection/ds_remove_failed/{pop}.ds_remove_failed.vcf.gz.tbi"
# 	benchmark:
# 		"benchmarks/{rule}/{pop}.{rule}.benchmark"
# 	params:
# 		prefix = "{pop}.{rule}"
# 	output:
# 		relatedness2 = "data/derived_data/sample_selection/{rule}/{pop}.{rule}.relatedness2",
# 		log = "data/derived_data/sample_selection/{rule}/{pop}.{rule}.log"
# 	shell:
# 		"/cluster/rcss-spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/vcftools-v0.1.14-ncnmzrsebn2ypj2ios5bfjzjmhwzooce/bin/vcftools --gzvcf {input.vcf} --out data/derived_data/sample_selection/{rule}/{params.prefix} --relatedness2"
