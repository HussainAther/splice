import pandas as pd

sampletable = pd.read_csv('config/sampletable.tsv', header=0, sep='\t')
sample_dict = pd.DataFrame.to_dict(sampletable) # Dictionary {"experiment" : [...], "control" : [...], "label" : [...]}
sample = sampletable["samplename"]
fullname = sampletable["orig_filename"]
control = sampletable.loc[sampletable["group"] == "control"]["orig_filename"]
treatment = sampletable.loc[sampletable["group"] == "treatment"]["orig_filename"]

def get_treatment_bams():
    return treatment

def get_control_bams():
    return control

def get_sample():
	return sample

rank = ["fdr_bottom5", "incleveldifference_top5", "incleveldifference_bottom5", "pvalue_bottom5"]

as_type = ["SE", "A5SS", "A3SS", "RI", "MXE"]

jc_type = ["JunctionCountOnly", "ReadsOnTargetAndJunctionCounts"]

sashimi = expand("rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
    rank=rank,
	sample = get_sample(),
    as_type=as_type,
    jc_type=jc_type)

bed12 = expand("rmats_out/{sample}/{as_type}.MATS.{jc_type}.bed",
        as_type=as_type,
		sample=get_sample(),
        jc_type=jc_type)

rule all:
    input: sashimi + bed12 

ruleorder: preparse > parse

rule rmats:
    input:
        gtf = "/data/athersh/Danio_rerio.GRCz11.94.gtf",
        b1 = get_treatment_bams(),
        b2 = get_control_bams(),
    output:
        expand("rmats_out/{sample}/summary.txt", sample = get_sample())
    params:
        length = 101,
        reading = "single"
    log:
        "logs/rmats"
    wrapper:
        "file:wrappers/rmats"

rule preparse:
	input: 
		"rmats_out/{sample}/summary.txt",
	output:
		"rmats_out/{sample}/{as_type}.MATS.{jc_type}.txt",
	wrapper:
		"file:wrappers/preparse"

rule parse:
    input:
        "rmats_out/{sample}/{as_type}.MATS.{jc_type}.txt",
    output:
        pvalue_bottom5 = "rmats_out/{sample}/pvalue_bottom5/{as_type}.MATS.{jc_type}.txt",
        fdr_bottom5 = "rmats_out/{sample}/fdr_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_bottom5 = "rmats_out/{sample}/incleveldifference_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_top5 = "rmats_out/{sample}/incleveldifference_top5/{as_type}.MATS.{jc_type}.txt",
    wrapper:
        "file:wrappers/parse"

rule sashimiplot:
    input: 
        rmats = expand("rmats_out/{sample}/{as_type}.MATS.{jc_type}.txt",sample=get_sample(), jc_type=jc_type, as_type=as_type),
        b1 = get_treatment_bams(),
        b2 = get_control_bams()
    output:
        event = expand("rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/Sashimi_index/{as_type}.event.list.txt", rank=rank, sample=get_sample(), jc_type=jc_type, as_type=as_type)
    params:
        dir = expand("rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/",sample=get_sample(), rank=rank, as_type=as_type, jc_type=jc_type),
        as_ =  expand("{as_type}",as_type=as_type)
#    log:
#        "logs/sashimiplot/{sample}/{rank}/{as_type}/{jc_type}"
    wrapper:
        "file:wrappers/sashimiplot"

rule pdfunite:
	input:
		"rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/Sashimi_index/{as_type}.event.list.txt",
	output: 
		"rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
	params:
		"rmats_out/{sample}/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/",
	wrapper:
		"file:wrappers/pdfunite"

rule convert_to_bed12:
    """
    For each output combination from rMATS, convert the events file to bed12 format
    """
    input: "rmats_out/{sample}/{as_type}.MATS.{jc_type}.txt"
    output: "rmats_out/{sample}/{as_type}.MATS.{jc_type}.bed"
    wrapper: "file:wrappers/convert_to_bed12"


# vim: ft=python
