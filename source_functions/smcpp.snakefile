#snakemake -s source_functions/smcpp.snakefile --jobs 6 --latency-wait 60 --rerun-incomplete --cluster-config source_functions/cluster/smcpp.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -n {cluster.n} --mem {cluster.mem} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/smcpp/190711.2.smcpp.log

import pandas as pd

configfile : "source_functions/config/smcpp.config.yaml"

# Dictionary of poplulation membership for all samples
samples_df = pd.read_csv(config['samples_file'])
# https://stackoverflow.com/questions/20112760/python-pandas-convert-dataframe-to-dictionary-with-multiple-values
samples_dict = {k: list(v) for k,v in samples_df.groupby("pop")["international_id"]}

# Need dictionary of randomly chosen distinguished samples
distinguished_df = pd.read_csv(config['distinguished_file'])
distinguished_dict = {k: list(v) for k,v in distinguished_df.groupby("pop")["distinguished"]}

rule all:
	input:
		expand("data/derived_data/smcpp/samples/{pop}_done.txt", pop = config['pop']),
		expand("data/derived_data/smcpp/plot_composite_estimate/plot_composite_estimate.{dataset}.{pop}.{rundate}.png",
		dataset = config['dataset'], pop = config['pop'], rundate = config['rundate'])

def string_choose(WC):
	string = WC.pop+":"+",".join(samples_dict[WC.pop])
	return string

# "This command will parse data for the contig chr1 for samples S1 and S2 which are members of population Pop1. You should run this once for each independent contig in your dataset, producing one SMC++ output file per contig."
rule vcf2smc:
	input:
		vcf = config['vcf'],
		index = config['vcf'] + ".tbi"
	params:
		chrom = "{chr}",
		# For the current population, string of all individuals in that population
		string = string_choose,
		# Select distinguished individual. Need to iterate over every chromosome of every possible distinguished individual
		distinguished = "{distinguished}",
		time_command = config['time_command'],
		time = 	"benchmarks/vcf2smc/vcf2smc.{dataset}.{pop}.{distinguished}.{chr}.{rundate}.time"
	output:
		smc = "data/derived_data/smcpp/vcf2smc/vcf2smc.{dataset}.{pop}.{distinguished}.{chr}.{rundate}.smc.gz"
	shell:
	# Specify input vcf, specify output directory, speficy chromosome, specify string of list of individuals in the population
	# "Note that "first" and "second" allele have no meaning for unphased data; if your data are not phased, it only makes sense to specify a single individual "
		"{params.time_command} {params.time} smc++ vcf2smc -d {params.distinguished} {params.distinguished} {input.vcf} {output.smc} {params.chrom} {params.string}"

# Generate population-level estimates
# 10/10/18 changed memory allocation in config file, uses a lot more memory than I previously thought
rule composite_estimate:
	input:
	#Need to use lambda function in order to correctly iterate over individuals in samples_dict
		 smc_chr = lambda wildcards: expand("data/derived_data/smcpp/vcf2smc/vcf2smc.{dataset}.{pop}.{distinguished}.{chr}.{rundate}.smc.gz", chr = config['chr'], dataset = wildcards.dataset, pop = wildcards.pop, distinguished = distinguished_dict[wildcards.pop], rundate = wildcards.rundate)
	params:
	# All chromosomes for all distinguished individuals
	# Should maybe sort into direcories like for estimate rules
		smc_in_string = "data/derived_data/smcpp/vcf2smc/vcf2smc.{dataset}.{pop}.*.{rundate}.smc.gz",
		model_out_dir = "data/derived_data/smcpp/composite_estimate/{dataset}.{pop}.{rundate}/",
		#thinning_k = thinning_choose,
		startdate = config['startdate'],
		nt = config['composite_estimate_nt'],
		mu = config['mu'],
		time_command = config['time_command'],
		time = "benchmarks/composite_estimate/composite_estimate.{dataset}.{pop}.{rundate}.time"
	output:
		model = "data/derived_data/smcpp/composite_estimate/{dataset}.{pop}.{rundate}/model.final.json"
	shell:
		"{params.time_command} {params.time} smc++ cv --cores {params.nt} --timepoints {params.startdate} 200000 -o {params.model_out_dir} {params.mu} {params.smc_in_string}"

#8/23/18 Jared asked me to run estimate for each distinguished individual and plot all individual lines together
#Potentially identify outlier "histories"
rule bootstrap_estimate:
	input:
	#Need to use lambda function in order to correctly iterate over individuals in samples_dict
		 smc_chr = lambda wildcards: expand('data/derived_data/smcpp/vcf2smc/vcf2smc.{dataset}.{pop}.{distinguished}.{chr}.{rundate}.smc.gz', dataset = wildcards.dataset, pop = wildcards.pop, distinguished = wildcards.distinguished, chr = config['chr'], rundate = wildcards.rundate)
	params:
		smc_in_dir = "data/derived_data/smcpp/vcf2smc/vcf2smc.{dataset}.{pop}.{distinguished}.*.{rundate}.smc.gz",
		model_out_dir = "data/derived_data/smcpp/bootstrap_estimate/{dataset}.{pop}.{rundate}/{distinguished}/",
		#thinning_k = thinning_choose,
		startdate = config['startdate'],
		# thinning_k = config['thinning_k'],
		nt = config['bootstrap_estimate_nt'],
		mu = config['mu'],
		time_command = config['time_command'],
		time = "benchmarks/bootstrap_estimate/bootstrap_estimate.{dataset}.{pop}.{distinguished}.{rundate}.time"
	output:
		model = "data/derived_data/smcpp/bootstrap_estimate/{dataset}.{pop}.{rundate}/{distinguished}/model.final.json"
	shell:
		"{params.time_command} {params.time} smc++ cv --cores 12 --timepoints {params.startdate} 200000 -o {params.model_out_dir} {params.mu} {params.smc_in_dir}"


rule plot_bootstrap_estimate:
	input:
		model = "data/derived_data/smcpp/bootstrap_estimate/{dataset}.{pop}.{rundate}/{distinguished}/model.final.json"
	output:
		plot = "data/derived_data/smcpp/plot_bootstrap_estimate/plot_bootstrap_estimate.{dataset}.{pop}.{distinguished}.{rundate}.png"
	shell:
		"smc++ plot {output.plot} {input.model} -g 5 -c"

rule plot_composite_estimate:
	input:
		model = "data/derived_data/smcpp/composite_estimate/{dataset}.{pop}.{rundate}/model.final.json"
	output:
		plot = "data/derived_data/smcpp/plot_composite_estimate/plot_composite_estimate.{dataset}.{pop}.{rundate}.png"
	#-k also plot spline knots
	shell:
		"smc++ plot {output.plot} {input.model} -g 5 -c"

rule done:
	input:
		bootstrap_plot = lambda wildcards: expand("data/derived_data/smcpp/plot_bootstrap_estimate/plot_bootstrap_estimate.{dataset}.{pop}.{distinguished}.{rundate}.png", dataset = config['dataset'], pop = wildcards.pop, distinguished = distinguished_dict[wildcards.pop], rundate = config['rundate'])
	output:
		done = "data/derived_data/smcpp/samples/{pop}_done.txt"
	shell:
		"echo 'done' > {output.done}"
