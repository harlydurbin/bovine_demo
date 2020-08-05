# snakemake -s source_functions/smartpca.snakefile -j 1000 --rerun-incomplete --keep-going --latency-wait 30 --use-conda --resources load=100 --config --cluster-config source_functions/cluster/smartpca.cluster.json --cluster "sbatch -p {cluster.p} -o {cluster.o} --account {cluster.account} -t {cluster.t} -c {cluster.c} --mem {cluster.mem} --account {cluster.account} --mail-user {cluster.mail-user} --mail-type {cluster.mail-type} --qos {cluster.qos}" -p &> log/snakemake_log/smartpca/200805.smartpca.log

import os
import re

include: "plink_qc.snakefile"

configfile: "source_functions/config/smartpca.config.yaml"

#Make log directories if they don't exist
os.makedirs("log/slurm_out/smartpca", exist_ok = True)
for x in expand("log/slurm_out/smartpca/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

os.makedirs("log/psrecord/smartpca", exist_ok = True)
for x in expand("log/psrecord/smartpca/{rules}", rules = config['rules']):
    os.makedirs(x, exist_ok = True)

rule smartpca_all:
	input:
		expand("data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.evec", thin_p = config['thin_p'], dataset = config['dataset'])

rule recode_input:
	input:
	# If projecting aurochs, need to use the use the same PLINK files but dataset name won't match, pull 'base_dataset' from dictionary in config
		bed = lambda wildcards: expand("data/derived_data/plink_qc/thin_variants/merge_thinned.{base_dataset}.{thin_p}.bed", base_dataset = config['dataset_key'][wildcards.dataset], thin_p = wildcards.thin_p),
		bim = lambda wildcards: expand("data/derived_data/plink_qc/thin_variants/merge_thinned.{base_dataset}.{thin_p}.bim", base_dataset = config['dataset_key'][wildcards.dataset], thin_p = wildcards.thin_p),
		fam = lambda wildcards: expand("data/derived_data/plink_qc/thin_variants/merge_thinned.{base_dataset}.{thin_p}.fam", base_dataset = config['dataset_key'][wildcards.dataset], thin_p = wildcards.thin_p)
	# If projecting aurochs, create a 'poplist' file without aurochs IDs
	params:
		poplist = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.poplist.txt"
	output:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind"
	run:
		if re.match('project', wildcards.dataset):
			shell('sed "s/-9/pheno/g" {input.fam} > {output.pedind}')
			shell('cp {input.bed} {output.bed}')
			shell('cp {input.bim} {output.pedsnp}')
			shell('cut -f1 {input.fam} | sed "/\\b\\(A2494_Yenisei_River\\|ancient_derbyshire\\)\\b/d" > {params.poplist}')
		else:
			shell('sed "s/-9/pheno/g" {input.fam} > {output.pedind}')
			shell('cp {input.bed} {output.bed}')
			shell('cp {input.bim} {output.pedsnp}')

rule create_par:
	input:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind"
	params:
		evec = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.evec",
		eval = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.eval",
		poplist = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.poplist.txt"
	output:
		par = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par"
	# Parameter file differs if aurochs are being projected onto PC space or included calculation
	run:
		if re.match('project', wildcards.dataset):
			shell('echo -e "genotypename:\\t{input.bed}\\nsnpname:\\t{input.pedsnp}\\nindivname:\\t{input.pedind}\\nevecoutname:\\t{params.evec}\\nevaloutname:\\t{params.eval}\\nnumchrom:\\t31\\nnumoutlieriter:\\t0\\nfastmode:\\tYES\\npoplist:\\t{params.poplist}" > {output.par}')
		else:
			shell('echo -e "genotypename:\\t{input.bed}\\nsnpname:\\t{input.pedsnp}\\nindivname:\\t{input.pedind}\\nevecoutname:\\t{params.evec}\\nevaloutname:\\t{params.eval}\\nnumchrom:\\t31\\nnumoutlieriter:\\t0\\nfastmode:\\tYES" > {output.par}')

rule smartpca:
	input:
		bed = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.bed",
		pedsnp = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedsnp",
		pedind = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.pedind",
		par = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.par"
	params:
		eigensoft_module = config['eigensoft_module'],
		psrecord = "log/psrecord/smartpca/smartpca/smartpca.{dataset}.{thin_p}.log"
	output:
		evec = "data/derived_data/smartpca/{dataset}.{thin_p}/smartpca.{dataset}.{thin_p}.evec"
	shell:
		"""
		module load {params.eigensoft_module}
		psrecord "smartpca -p {input.par}" --log {params.psrecord} --include-children --interval 5
		"""
