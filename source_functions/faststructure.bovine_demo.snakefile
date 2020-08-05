# snakemake -s source_functions/faststructure.bovine_demo.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --resources load=100 --config --cluster-config source_functions/cluster/faststructure.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/faststructure/200805.faststructure.log

include: "plink_qc.snakefile"

configfile: "source_functions/config/faststructure.config.yaml"

# Make log directories if they don't exist
os.makedirs("log/slurm_out/faststructure", exist_ok = True)
for x in expand("log/slurm_out/faststructure/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/faststructure", exist_ok = True)
for x in expand("log/psrecord/faststructure/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

rule faststructure_all:
	input:
		expand("data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}.ML.txt", thin_p = config['thin_p'], dataset = config['dataset'])

rule structure:
	input:
		merged_bed = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.bed",
		merged_bim = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.bim",
		merged_fam = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}.fam"
	# Snakemake looks for envs/ relative to the snakefile
	conda:
		"envs/faststructure.yaml"
	params:
		in_prefix = "data/derived_data/plink_qc/thin_variants/merge_thinned.{dataset}.{thin_p}",
		out_prefix = "data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}",
		k = "{k}"
	output:
		structure_files = expand("data/derived_data/faststructure/structure/{{dataset}}.{{thin_p}}/structure.{{dataset}}.{{thin_p}}.{{k}}.{extension}", extension = ["meanP", "meanQ", "varP", "varQ"])
	shell:
		"""
		structure.py -K {params.k} --input={params.in_prefix} --output={params.out_prefix} --prior=simple --cv=0 --full
		"""

rule structure_ml:
	input:
		structure_files = lambda wildcards: expand("data/derived_data/faststructure/structure/{dataset}.{thin_p}/structure.{dataset}.{thin_p}.{k}.{extension}", dataset = wildcards.dataset, thin_p = wildcards.thin_p, k = config['k'][wildcards.dataset], extension = ["meanP", "meanQ", "varP", "varQ"])
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
