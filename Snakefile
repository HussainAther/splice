import pandas as pd

sampletable = pd.read_csv('config/sampletable.tsv', header=0, sep='\t')
sample_dict = pd.DataFrame.to_dict(sampletable) # Dictionary {"experiment" : [...], "control" : [...], "label" : [...]}
sample = sampletable["samplename"]
fullname = sampletable["orig_filename"]
control = sampletable.loc[sampletable["group"] == "control"]["samplename"]
treatment = sampletable.loc[sampletable["group"] == "treatment"]["samplename"]

def get_treatment_bams():
    return treatment

def get_control_bams():
    return control

scientist = ["DC_shep_kd", "JM_rump_kd", "JM_rump_overexpression"]

rank = ["fdr_bottom5", "incleveldifference_top5", "incleveldifference_bottom5", "pvalue_bottom5"]

as_type = ["SE", "A5SS", "A3SS", "RI", "MXE"]

jc_type = ["JunctionCountOnly", "ReadsOnTargetAndJunctionCounts"]

config["gtf"] = "/data/Lei_student/Hussain/RNASeq/features.gtf"

outdir = "rmats_out"

sashimi = expand("{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
    outdir=outdir,
    scientist=scientist,
    rank=rank,
    as_type=as_type,
    jc_type=jc_type)

bed12 = expand("{outdir}/{scientist}/MATS_output/{as_type}.MATS.{jc_type}.bed",
        outdir=outdir,
        scientist=scientist,
        as_type=as_type,
        jc_type=jc_type)

subworkflow rnaseq:
    """
    Retrieve files from the RNA-Seq workflow
    """
    workdir: "../lcdb-wf/lcdb-wf/workflows/rnaseq/"
    snakefile: "../lcdb-wf/lcdb-wf/workflows/rnaseq/Snakefile"

bigwig = rnaseq(expand('data/rnaseq_samples/{sample}/{sample}.cutadapt.bam.pos.bigwig', sample=sample)) + rnaseq(expand('data/rnaseq_samples/{sample}/{sample}.cutadapt.bam.neg.bigwig', sample=sample))

bigwig = expand("{sample}.neg.bigwig", sample=fullname) + expand("{sample}.pos.bigwig", sample=fullname)

rule all:
    input: sashimi + bed12 

rule rmats:
    input:
        gtf = config["gtf"],
        b1 = get_treatment_bams(),
        b2 = get_control_bams()
    output:
        "{outdir}/{scientist}/MATS_output/summary.txt"
    params:
        length = 51,
        reading = "single"
    log:
        "logs/rmats"
    wrapper:
        "file:wrappers/rmats"

rule parse:
    input:
        "{outdir}/{scientist}/MATS_output/{as_type}.MATS.{jc_type}.txt",
    output:
        pvalue_bottom5 = "{outdir}/{scientist}/MATS_output/pvalue_bottom5/{as_type}.MATS.{jc_type}.txt",
        fdr_bottom5 = "{outdir}/{scientist}/MATS_output/fdr_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_bottom5 = "{outdir}/{scientist}/MATS_output/incleveldifference_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_top5 = "{outdir}/{scientist}/MATS_output/incleveldifference_top5/{as_type}.MATS.{jc_type}.txt",
    wrapper:
        "file:wrappers/parse"

rule sashimiplot:
    input: 
        rmats = "{outdir}/{scientist}/MATS_output/{rank}/{as_type}.MATS.{jc_type}.txt",
        b1 = get_treatment_bams(),
        b2 = get_control_bams()
    output:
        dir = "{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/",
        subdir = "{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/",
        event = "{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_index/{as_type}.event.list.txt"
    params:
        as_ =  "{as_type}",
    log:
        "logs/sashimiplot"
    wrapper:
        "file:wrappers/sashimiplot"

rule pdfunite:
    input: 
        "{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/",
    output: 
        "{outdir}/{scientist}/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
    log:
        "logs/pdfunite"
    shell:
        "pdfunite {input}* {output} &> {log}"

rule convert_to_bed12:
    """
    For each output combination from rMATS, convert the events file to bed12 format
    """
    input: "{outdir}/{scientist}/MATS_output/{as_type}.MATS.{jc_type}.txt"
    output: "{outdir}/{scientist}/MATS_output/{as_type}.MATS.{jc_type}.bed"
    wrapper: "file:wrappers/convert_to_bed12"


# vim: ft=python
