# snakemake -s source_functions/plink_qc.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --config --cluster-config source_functions/cluster/plink_qc.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type}" -p &> log/snakemake_log/joint_genotyping/200804.plink_qc.log

import os

# include path is relative to the path of this file
include: "post_process.snakefile"

include: "sex_imputation.snakefile"

configfile: "source_functions/config/plink_qc.yaml"

os.makedirs("log/slurm_out/plink_qc", exist_ok = True)

# Make log directories if they don't exist
for x in expand("log/slurm_out/plink_qc/{rules}", rules = config['plink_qc_rules']):
    os.makedirs(x, exist_ok = True)

for x in expand("log/psrecord/joint_genotyping/{rules}", rules = config['plink_qc_rules']):
    os.makedirs(x, exist_ok = True)

rule plink_qc_all:
	input:
		expand("data/derived_data/plink_qc/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.{extension}", downsample_dataset = config['downsample_dataset'], thin_p = config['thin_p'], extension = ['bed', 'bim', 'fam']), expand("data/derived_data/plink_qc/thin_variants/merge_thinned.full.{thin_p}.{extension}", thin_p = config['thin_p'], extension = ['bed', 'bim', 'fam'])

rule qc_autosomes:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.{autosome}.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		prefix = "data/derived_data/plink_qc/initial_qc/initial_qc.{autosome}",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter'],
		autosome = "{autosome}"
	output:
		bed = "data/derived_data/plink_qc/initial_qc/initial_qc.{autosome}.bed",
		bim = "data/derived_data/plink_qc/initial_qc/initial_qc.{autosome}.bim",
		fam = "data/derived_data/plink_qc/initial_qc/initial_qc.{autosome}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr {params.autosome} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

rule qc_x:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.X.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		pab = config['pab'],
		prefix = "data/derived_data/plink_qc/initial_qc/initial_qc.X",
		nt = config['plink_nt'],
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/plink_qc/initial_qc/initial_qc.X.bed",
		bim = "data/derived_data/plink_qc/initial_qc/initial_qc.X.bim",
		fam = "data/derived_data/plink_qc/initial_qc/initial_qc.X.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr X -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --from-bp 1 --to-bp {params.pab} --update-sex {input.imputed_sexes} --geno {params.geno_filter} --out {params.prefix}
		"""

# set heterozygous calls to missing on the Y
rule qc_y:
	input:
		vcf = "data/derived_data/joint_genotyping/remove_samples/bovine_demo.Y.vcf.gz",
		imputed_sexes = config['imputed_sexes']
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		prefix = "data/derived_data/plink_qc/initial_qc/initial_qc.Y",
		geno_filter = config['geno_filter']
	output:
		bed = "data/derived_data/plink_qc/initial_qc/initial_qc.Y.bed",
		bim = "data/derived_data/plink_qc/initial_qc/initial_qc.Y.bim",
		fam = "data/derived_data/plink_qc/initial_qc/initial_qc.Y.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --vcf {input.vcf} --make-bed --double-id --cow --chr Y --set-hh-missing --geno {params.geno_filter} -set-missing-var-ids @:#\$1\$2 --threads {params.nt} --update-sex {input.imputed_sexes} --out {params.prefix}
		"""

def thin_chooser(WC):
	thinness = config["thin_p_dict"][WC.thin_p] # References the dictionary, then the value associated with that key (which is your wildcard)
	return thinness # Returns this, which is a param in your rule, which is inturn pluged into the shell command

# Output is all samples, variant density thinned to specified thinning parameter, chromosome-by-chromosome
# Can't use update-sex and update-ids in the same run, do it here instead
rule thin_variants:
	input:
		bed = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.bed",
		bim = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.bim",
		fam = "data/derived_data/plink_qc/initial_qc/initial_qc.{chr}.fam",
		rename_file = config['rename_file']
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
		plink --bfile {params.in_prefix} --maf .0001 --thin {params.thin_p} --update-ids {input.rename_file} --make-bed --double-id --cow --threads {params.nt} --out {params.out_prefix}
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
