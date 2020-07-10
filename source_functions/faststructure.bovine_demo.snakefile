# snakemake -s source_functions/madtom.faststructure.snakefile -j 4 --rerun-incomplete --latency-wait 30 --config -p &> log/snakemake_log/200325.madtom.faststructure.log

configfile: "source_functions/config/faststructure.config.yaml"


#Make log directories if they don't exist
os.makedirs("log/slurm_out/faststructure", exist_ok = True)
for x in expand("log/slurm_out/faststructure/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/faststructure", exist_ok = True)

rule all:
	input:
		expand("data/derived_data/faststructure/structure/{dataset}_{thin_p}/structure.{dataset}_{thin_p}.ML.txt", thin_p = config['thin_p'], dataset = config['downsample_dataset'] + ['all'])

def thin_chooser(WC):
	thinness = config["thin_p_dict"][WC.thin_p] # References the dictionary, then the value associated with that key (which is your wildcard)
	return thinness # Returns this, which is a param in your rule, which is inturn pluged into the shell command

rule thin_variants:
	input:
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{chr}.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{chr}.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{chr}.fam"
	params:
		plink_module = config['plink_module'],
		plink_nt = config['plink_nt'],
		in_prefix = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.{chr}",
		out_prefix = "data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{thin_p}",
		thin_p = thin_chooser
	output:
		thinned_bed = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{thin_p}.bed"),
		thinned_bim = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{thin_p}.bim"),
		thinned_fam = temp("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{thin_p}.fam")
	shell:
		"""
		module load {param.plink_module}
		plink --bfile {params.in_prefix} --mac 1 --thin {params.thin_p} --make-bed --double-id --cow --threads {params.nt} --out {params.out_prefix}
		"""

rule merge_thinned:
	input:
		thinned_bed = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{{thin_p}}.bed", chr = config['chr']),
		thinned_bim = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{{thin_p}}.bim", chr = config['chr']),
		thinned_fam = expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{{thin_p}}.fam", chr = config['chr'])
	params:
		plink_module = config['plink_module'],
		plink_nt = config['plink_nt'],
		prefixes = lambda wildcards: expand("data/derived_data/faststructure/thin_variants/thin_variants.{chr}_{thin_p}\n", thin_p = wildcards.thin_p, chr = config['chr']),
		out_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}"
	output:
		merge_list = "data/derived_data/faststructure/thin_variants/merge_list.all_{thin_p}.list",
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.fam"
	shell:
		"""
		module load {param.plink_module}
		echo -e '{params.prefixes}' | sed 's/^ *//g' > {output.merge_list}
		plink --merge_list {output.merge_list} --double-id --cow --threads {params.nt} --make-bed --out {params.out_prefix}
		"""

rule downsample_indiv:
	input:
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.all.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}.fam",
		keep_list = "data/derived_data/faststructure/thin_variants/keeplist.{downsample_dataset}.txt"
	params:
		plink_module = config['plink_module'],
		plink_nt = config['plink_nt'],
		in_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.all_{thin_p}",
		out_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}_{thin_p}"
	output:
		downsample_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}_{thin_p}.bed",
		downsample_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}_{thin_p}.bim",
		downsample_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.{downsample_dataset}_{thin_p}.fam"
	shell:
		"""
		module load {param.plink_module}
		plink --bfile {params.in_prefix} --double-id --cow --threads {params.nt} --keep {params.keep_list} --make-bed --out {params.out_prefix}
		"""

rule structure:
	input:
		merged_bed = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}_{thin_p}.bed",
		merged_bim = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}_{thin_p}.bim",
		merged_fam = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}_{thin_p}.fam"
	conda:
		"source_functions/envs/faststructure.yaml"
	params:
		in_prefix = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}_{thin_p}",
		out_prefix = "data/derived_data/faststructure/structure/{dataset}_{thin_p}/structure.{dataset}_{thin_p}.k{k}",
		k = "{k}",
		psrecord = "log/psrecord/faststructure/structure/structure.{dataset}_{thin_p}.k{k}.log"
	output:
		structure_files = expand("data/derived_data/faststructure/structure/{{dataset}}_{{thin_p}}/structure.{{dataset}}_{{thin_p}}.k{{k}}.{suffix}", suffix = ["meanP", "meanQ", "varP", "varQ"])
	shell:
		"""
		structure.py -K {params.k} --input={params.prefix} --output={params.prefix} --prior=simple --cv=0 --full
		"""


rule structure_ml:
	input:
		structure_files = expand("data/derived_data/faststructure/structure/{{dataset}}_{{thin_p}}/structure.{{dataset}}_{{thin_p}}.k{k}.{suffix}", k = config['k'], suffix = ["meanP", "meanQ", "varP", "varQ"])
	conda:
		"source_functions/envs/faststructure.yaml"
	params:
		prefix = "data/derived_data/faststructure/structure/{dataset}_{thin_p}/structure.{dataset}_{thin_p}"
	output:
		ML = "data/derived_data/faststructure/structure/{dataset}_{thin_p}/structure.{dataset}_{thin_p}.ML.txt"
	shell:
		"""
		chooseK.py --input={params.prefix} &> {output.ML}
		"""
