# snakemake -s source_functions/smartpca.bovine_demo.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --resources load=100 --config --cluster-config source_functions/cluster/smartpca.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/smartpca/200713.smartpca.log

configfile: "source_functions/config/smartpca.config.yaml"

#Make log directories if they don't exist
os.makedirs("log/slurm_out/smartpca", exist_ok = True)
for x in expand("log/slurm_out/smartpca/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/smartpca", exist_ok = True)

rule all:
	input:
		expand("data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par", thin_p = config['thin_p'], dataset = config['dataset'])

# fam = pedind
# bim = pedsnp
# bed = bed

rule recode_input:
	input:
		# bed = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.bed",
		# bim = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.bim",
		# fam = "data/derived_data/faststructure/thin_variants/merge_thinned.{dataset}.{thin_p}.fam"
		bed = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.1.bed",
		bim = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.1.bim",
		fam = "data/derived_data/joint_genotyping/phasing_qc/phasing_qc.1.fam"
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
