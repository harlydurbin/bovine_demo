# snakemake -s source_functions/treemix.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --config --cluster-config source_functions/cluster/treemix.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/treemix/200817.treemix.log

configfile: "source_functions/config/treemix.yaml"

# Make log directories if they don't exist
os.makedirs("log/slurm_out/treemix", exist_ok = True)
for x in expand("log/slurm_out/treemix/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/treemix", exist_ok = True)
os.makedirs("log/psrecord/treemix/treemix", exist_ok = True)

rule target:
	input:
		targ = lambda wildcards: expand("data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.vertices.gz", dataset = config['dataset'], thin_p = config['thin_p'], m = list(range(config['min_m'], config['max_m'])))

rule define_cluster:
	input:
		fam = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.fam",
		script = "source_functions/define_cluster.R"
	params:
		r_module = config['r_module']
	output:
		within = "data/derived_data/treemix/define_cluster.{dataset}.{thin_p}.txt"
	shell:
		"""
		module load {params.r_module}
		Rscript --vanilla {input.script} {input.fam} {output.within}
		"""

rule write_cluster:
	input:
		bed = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bed",
		bim = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bim",
		fam = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.fam",
		within = "data/derived_data/treemix/define_cluster.{dataset}.{thin_p}.txt"
	params:
		prefix = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}",
		plink_module = config['plink_module'],
	output:
		clst = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.clst"
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.prefix} --within {input.within} --double-id --cow --write-cluster --out {params.prefix}
		"""

rule freq:
	input:
		bed = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bed",
		bim = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bim",
		fam = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.fam",
		clst = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.clst"
	params:
		in_prefix = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}",
		out_prefix = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}",
		plink_module = config['plink_module'],
		strat = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.strat"
	output:
		strat_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.strat.gz"
	shell:
		"""
		module load {params.plink_module}
		plink --bfile {params.in_prefix} --freq --double-id --cow --within {input.clst} --out {params.out_prefix}
		gzip {params.strat}
		"""

# https://anaconda.org/bioconda/treemix
rule plink2tm:
	input:
		strat_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.strat.gz",
		script = "source_functions/plink2treemix.py"
	# script is written in python2.7
	params:
		python_module = config['python_module']
	output:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.gz"
	shell:
		"""
		module load {params.python_module}
		{input.script} {input.strat_gz} {output.treemix_gz}
		"""

rule treemix_base:
	input:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.gz"
	params:
		ld_snps = config['ld_snps'],
		root = config['root'],
		prefix = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base",
		psrecord = "log/psrecord/treemix/treemix/treemix.{dataset}.{thin_p}.base.log"
	output:
		cov_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.cov.gz",
		covse_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.covse.gz",
		modelcov_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.modelcov.gz",
		treeout_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.treeout.gz",
		vertices_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.vertices.gz",
		edges_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.edges.gz"
	shell:
		"""
		psrecord "treemix -i {input.treemix_gz} -k {params.ld_snps} -se -root {params.root} -o {params.prefix}" --log {params.psrecord} --include-children --interval 5
		"""

rule treemix:
	input:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.gz",
		cov_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.cov.gz",
		covse_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.covse.gz",
		modelcov_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.modelcov.gz",
		treeout_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.treeout.gz",
		vertices_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.vertices.gz",
		edges_in = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.base.edges.gz"
	params:
		ld_snps = config['ld_snps'],
		root = config['root'],
		prefix = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}",
		m = "{m}",
		psrecord = "log/psrecord/treemix/treemix/treemix.{dataset}.{thin_p}.{m}.log"
	output:
		cov_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.cov.gz",
		covse_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.covse.gz",
		modelcov_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.modelcov.gz",
		treeout_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.treeout.gz",
		vertices_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.vertices.gz",
		edges_out = "data/derived_data/treemix/output/{dataset}.{thin_p}/treemix.{dataset}.{thin_p}.{m}.edges.gz"
	shell:
		"""
		psrecord "treemix -i {input.treemix_gz} -k {params.ld_snps} -se -root {params.root} -m {params.m} -g {input.vertices_in} {input.edges_in} -o {params.prefix}" --log {params.psrecord} --include-children --interval 5
		"""
