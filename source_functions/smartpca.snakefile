# snakemake -s source_functions/smartpca.bovine_demo.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --resources load=100 --config --cluster-config source_functions/cluster/smartpca.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/smartpca/200713.smartpca.log

include: "plink_qc.snakefile"

configfile: "source_functions/config/smartpca.config.yaml"

#Make log directories if they don't exist
os.makedirs("log/slurm_out/smartpca", exist_ok = True)
for x in expand("log/slurm_out/smartpca/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/slurm_out/plink_qc/{qc_rules}", qc_rules = ['thin_variants', 'merge_thinned', 'downsample_indiv']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/smartpca", exist_ok = True)

rule smartpca_all:
	input:
		expand("data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par", thin_p = config['thin_p'], dataset = config['dataset'])

def thin_chooser(WC):
	thinness = config["thin_p_dict"][WC.thin_p] # References the dictionary, then the value associated with that key (which is your wildcard)
	return thinness # Returns this, which is a param in your rule, which is inturn pluged into the shell command

# Output is all samples, variant density thinned to specified thinning parameter, chromosome-by-chromosome
rule thin_variants:
	input:
		bed = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.bed",
		bim = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.bim",
		fam = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.fam"
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		in_prefix = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}",
		out_prefix = "data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{thin_p}",
		thin_p = thin_chooser
	resources:
		load = 1
	output:
		thinned_bed = "data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{thin_p}.bed",
		thinned_bim = "data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{thin_p}.bim",
		thinned_fam = "data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{thin_p}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.in_prefix} --maf .0001 --thin {params.thin_p} --make-bed --double-id --cow --threads {params.nt} --out {params.out_prefix}
		"""

# Input is all chromosomes for the specified thinning parameter
# Output is all samples, variant density thinned to specified thinning parameter, merged to whole-genome
rule merge_thinned:
	input:
		thinned_bed = expand("data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{{thin_p}}.bed", chr = config['chr']),
		thinned_bim = expand("data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{{thin_p}}.bim", chr = config['chr']),
		thinned_fam = expand("data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{{thin_p}}.fam", chr = config['chr'])
	resources:
		load = 1
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		prefixes = lambda wildcards: expand("data/derived_data/plink_qc/thin_variants/thin_variants.{chr}.full.{thin_p}\n", thin_p = wildcards.thin_p, chr = config['chr']),
		out_prefix = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}"
	output:
		merge_list = "data/derived_data/plink_qc/thin_variants/merge_list.full.{thin_p}.list",
		merged_bed = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.bed",
		merged_bim = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.bim",
		merged_fam = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.fam"
	shell:
		"""
		module load {params.plink_module}
		echo -e '{params.prefixes}' | sed 's/^ *//g' > {output.merge_list}
		plink --merge-list {output.merge_list} --double-id --cow --threads {params.nt} --make-bed --out {params.out_prefix}
		"""

# Output is all samples, variant density thinned to specified thinning parameter, merged to whole-genome, only samples in specified downsample_dataset
rule downsample_indiv:
	input:
		merged_bed = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.bed",
		merged_bim = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.bim",
		merged_fam = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.fam",
		keep_list = lambda wildcards: expand("data/derived_data/plink_qc/thin_variants/keeplist.{downsample_dataset}.txt", downsample_dataset = wildcards.downsample_dataset)
	resources:
		load = 1
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		in_prefix = "data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}",
		out_prefix = "data/derived_data/plink_qc/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}"
	output:
		downsample_bed = "data/derived_data/plink_qc/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.bed",
		downsample_bim = "data/derived_data/plink_qc/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.bim",
		downsample_fam = "data/derived_data/plink_qc/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.in_prefix} --double-id --cow --threads {params.nt} --keep {input.keep_list} --maf .0001 --make-bed --out {params.out_prefix}
		"""

# fam = pedind
# bim = pedsnp
# bed = bed

rule recode_input:
	input:
		bed = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.bed",
		bim = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.bim",
		fam = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.fam"
	output:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind"
	shell:
		"""
		sed 's/-9/pheno/g' {input.fam} > {output.pedind}
		cp {input.bed} {output.bed}
		cp {input.bim} {output.pedsnp}
		"""

rule create_par:
	input:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind"
	params:
		evec = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.evec",
		eval = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.eval"
	output:
		par = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par"
	shell:
		"""
		echo -e "genotypename:\\t{input.bed}\\nsnpname:\\t{input.pedsnp}\\nindivname:\\t{input.pedind}\\nevecoutname:\\t{params.evec}\\nevaloutname:\\t{params.eval}\\nnumchrom:\\t31\\nfastmode:\\tYES" > {output.par}
		"""

rule smartpca:
	input:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind",
		par = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par"
	params:
		eigensoft_module = config['eigensoft_module'],
		psrecord = psrecord = "log/psrecord/smartpca/smartpca/smartpca.{dataset}.{thin_p}.log"
	output:
		evec = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.evec",
		eval = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.eval"
	shell:
		"""
		module load {params.eigensoft_module}
		psrecord "smartpca -p {input.par}" --log {params.psrecord} --include-children --interval 5
		"""
