#/CIFS/MUG01_N/deckerje/tnr343/dominica/filtered/merged.thinned_dominica.*

#snakemake -s source_functions/treemix.snakefile --rerun-incomplete --latency-wait 120 --jobs 10 &> snakemake_logs/treemix/191206.treemix.log

configfile: "source_functions/config/191204.treemix.json"

m_list = list(range(config['min_m'], config['max_m']))

rule target:
	input:
		targ = lambda wildcards: expand("data/derived_data/treemix/output/{dataset}/treemix.{dataset}.{m}.vertices.gz", dataset = config['dataset'], m = m_list)

rule keep:
	input:
		keepfile = "data/plink/{dataset}.keep.txt"
	params:
		prefix = "data/plink/{dataset}"
	output:
		bed = "data/plink/{dataset}.bed",
		bim = "data/plink/{dataset}.bim",
		fam = "data/plink/{dataset}.fam"
	shell:
		"plink --bfile data/plink/full --keep {input.keepfile} --make-bed --out {params.prefix}"

rule clst:
	input:
		bed = "data/plink/{dataset}.bed",
		bim = "data/plink/{dataset}.bim",
		fam = "data/plink/{dataset}.fam"
	params:
		prefix = "data/plink/{dataset}"
	output:
		clst = "data/plink/{dataset}.clst"
	shell: "plink --bfile {params.prefix} --write-cluster --family --out {params.prefix}"

rule freq:
	input:
		bed = "data/plink/{dataset}.bed",
		bim = "data/plink/{dataset}.bim",
		fam = "data/plink/{dataset}.fam",
		clst = "data/plink/{dataset}.clst"
	params:
		prefix_in = "data/plink/{dataset}",
		prefix_out = "data/derived_data/treemix/plink2tm/{dataset}"
	output:
		strat = "data/derived_data/treemix/plink2tm/{dataset}.frq.strat"
	shell: "plink --bfile {params.prefix_in} --freq --within {input.clst} --out {params.prefix_out}"

#plink --bfile data/plink/treemix_pops --freq --within data/plink/treemix_pops.clst --out data/derived_data/treemix/plink2tm/treemix_pops

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
