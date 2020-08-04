# snakemake -s source_functions/faststructure.bovine_demo.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --resources load=100 --config --cluster-config source_functions/cluster/faststructure.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/faststructure/200713.faststructure.log

include: "plink_qc.snakefile"

configfile: "source_functions/config/faststructure.config.yaml"

#Make log directories if they don't exist
os.makedirs("log/slurm_out/faststructure", exist_ok = True)
for x in expand("log/slurm_out/faststructure/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/faststructure", exist_ok = True)

rule faststructure_all:
	input:
		expand("data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}.ML.txt", thin_p = config['thin_p'], dataset = config['dataset'])

def thin_chooser(WC):
	thinness = config["thin_p_dict"][WC.thin_p] # References the dictionary, then the value associated with that key (which is your wildcard)
	return thinness # Returns this, which is a param in your rule, which is inturn pluged into the shell command

# Output is all samples, variant density thinned to specified thinning parameter, chromosome-by-chromosome
rule thin_variants:
	input:
		bed = "data/derived_data/joint_genotyping/plink_qc/qc.{chr}.bed",
		bim = "data/derived_data/joint_genotyping/plink_qc/qc.{chr}.bim",
		fam = "data/derived_data/joint_genotyping/plink_qc/qc.{chr}.fam"
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		in_prefix = "data/derived_data/joint_genotyping/plink_qc/qc.{chr}",
		out_prefix = "data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{thin_p}",
		thin_p = thin_chooser
	resources:
		load = 1
	output:
		thinned_bed = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{thin_p}.bed"),
		thinned_bim = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{thin_p}.bim"),
		thinned_fam = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{thin_p}.fam")
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.in_prefix} --maf .0001 --thin {params.thin_p} --make-bed --double-id --cow --threads {params.nt} --out {params.out_prefix}
		"""

# Input is all chromosomes for the specified thinning parameter
# Output is all samples, variant density thinned to specified thinning parameter, merged to whole-genome
rule merge_thinned:
	input:
		thinned_bed = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{{thin_p}}.bed", chr = config['chr']),
		thinned_bim = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{{thin_p}}.bim", chr = config['chr']),
		thinned_fam = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{{thin_p}}.fam", chr = config['chr'])
	resources:
		load = 1
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		prefixes = lambda wildcards: expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}.full.{thin_p}\n", thin_p = wildcards.thin_p, chr = config['chr']),
		out_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}"
	output:
		merge_list = "data/derived_data/faststructure/thin_variants/merge_list.full.{thin_p}.list",
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.fam"
	shell:
		"""
		module load {params.plink_module}
		echo -e '{params.prefixes}' | sed 's/^ *//g' > {output.merge_list}
		plink --merge_list {output.merge_list} --double-id --cow --threads {params.nt} --make-bed --out {params.out_prefix}
		"""

# Output is all samples, variant density thinned to specified thinning parameter, merged to whole-genome, only samples in specified downsample_dataset
rule downsample_indiv:
	input:
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}.fam",
		keep_list = "data/derived_data/faststructure/thin_variants/keeplist.{downsample_dataset}.txt"
	resources:
		load = 1
	params:
		plink_module = config['plink_module'],
		nt = config['plink_nt'],
		in_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.full.{thin_p}",
		out_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}"
	output:
		downsample_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.bed",
		downsample_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.bim",
		downsample_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}.{thin_p}.fam"
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.in_prefix} --double-id --cow --threads {params.nt} --keep {input.keep_list} --maf .0001 --make-bed --out {params.out_prefix}
		"""

rule structure:
	input:
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.fam"
	# Snakemake looks for envs/ relative to the snakefile
	conda:
		"envs/faststructure.yaml"
	resources:
		load = 10
	params:
		in_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}",
		out_prefix = "data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}.k{k}",
		k = "{k}",
		psrecord = "log/psrecord/faststructure/structure/structure.{dataset}.{thin_p}.k{k}.log"
	output:
		structure_files = expand("data/derived_data/faststructure/structure/{{dataset}}.{{thin_p}}/structure.{{dataset}}.{{thin_p}}.k{{k}}.{suffix}", suffix = ["meanP", "meanQ", "varP", "varQ"])
	shell:
		"""
		psrecord "structure.py -K {params.k} --input={params.in_prefix} --output={params.out_prefix} --prior=simple --cv=0 --full" --log {params.psrecord} --include-children --interval 5
		"""


rule structure_ml:
	input:
		structure_files = expand("data/derived_data/faststructure/structure/{{dataset}}.{{thin_p}}/structure.{{dataset}}.{{thin_p}}.k{k}.{suffix}", k = config['k'], suffix = ["meanP", "meanQ", "varP", "varQ"])
	resources:
		load = 1
	conda:
		"envs/faststructure.yaml"
	params:
		prefix = "data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}"
	output:
		ML = "data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}.ML.txt"
	shell:
		"""
		chooseK.py --input={params.prefix} &> {output.ML}
		"""
