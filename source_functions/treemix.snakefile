#/CIFS/MUG01_N/deckerje/tnr343/dominica/filtered/merged.thinned_dominica.*

#snakemake -s source_functions/treemix.snakefile --rerun-incomplete --latency-wait 120 --jobs 10 &> snakemake_logs/treemix/191206.treemix.log

configfile: "source_functions/config/treemix.json"

m_list = list(range(config['min_m'], config['max_m']))

rule target:
	input:
		targ = lambda wildcards: expand("data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.vertices.gz", dataset = config['dataset'], m = m_list)

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
		fam = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.fam"
	params:
		prefix = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}"
	output:
		clst = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.clst"
	shell: "plink --bfile {params.prefix} --within --write-cluster --out {params.prefix}"

rule freq:
	input:
		bed = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bed",
		bim = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.bim",
		fam = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.fam"
		clst = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}.clst"
	params:
		in_prefix = "data/derived_data/plink_qc/thin_variants/{dataset}.{thin_p}/merge_thinned.{dataset}.{thin_p}",
		out_prefix = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}"
	output:
		strat = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.{thin_p}.frq.strat"
	shell: "plink --bfile {params.in_prefix} --freq --within {input.clst} --out {params.out_prefix}"

rule gz:
	input:
		strat = "data/derived_data/treemix/plink2tm/{dataset}.frq.strat"
	output:
		strat_gz = "data/derived_data/treemix/plink2tm/{dataset}.frq.strat.gz"
	shell: "pigz {input.strat}"

rule plink2tm:
	input:
		strat_gz = "data/derived_data/treemix/plink2tm/{dataset}.frq.strat.gz",
		script = "source_functions/plink2treemix.py"
	output:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.frq.gz"
	shell: "python2.7 source_functions/plink2treemix.py {input.strat_gz} {output.treemix_gz}"

rule treemix_base:
	input:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.frq.gz"
	benchmark:
		"benchmarks/treemix/treemix.{dataset}.base.benchmark"
	params:
		ld_snps = config['ld_snps'],
		root = config['root'],
		prefix = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base"
	output:
		cov_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.cov.gz",
		covse_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.covse.gz",
		modelcov_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.modelcov.gz",
		treeout_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.treeout.gz",
		vertices_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.vertices.gz",
		edges_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.edges.gz"
	shell: "treemix -i {input.treemix_gz} -k {params.ld_snps} -se -root {params.root} -o {params.prefix}"

# def previous_tree(WC):
# 	oldm = str(int(WC.m) - 1)
# 	filestring = "data/derived_data/treemix/output/" + WC.dataset + "/treemix." + WC.dataset + "." + oldm + ".vertices.gz"
# 	return filestring

rule treemix:
	input:
		treemix_gz = "data/derived_data/treemix/plink2tm/plink2tm.{dataset}.frq.gz",
		cov_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.cov.gz",
		covse_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.covse.gz",
		modelcov_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.modelcov.gz",
		treeout_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.treeout.gz",
		vertices_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.vertices.gz",
		edges_in = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.base.edges.gz"
	benchmark:
		"benchmarks/treemix/treemix.{dataset}.{m}.benchmark"
	params:
		ld_snps = config['ld_snps'],
		root = config['root'],
		prefix = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}",
		m = "{m}"
	output:
		cov_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.cov.gz",
		covse_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.covse.gz",
		modelcov_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.modelcov.gz",
		treeout_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.treeout.gz",
		vertices_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.vertices.gz",
		edges_out = "data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.edges.gz"
	shell: "treemix -i {input.treemix_gz} -k {params.ld_snps} -se -root {params.root} -m {params.m} -g {input.vertices_in} {input.edges_in} -o {params.prefix}"
