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

rank = ["fdr_bottom5", "incleveldifference_top5", "incleveldifference_bottom5", "pvalue_bottom5"]

as_type = ["SE", "A5SS", "A3SS", "RI", "MXE"]

jc_type = ["JunctionCountOnly", "ReadsOnTargetAndJunctionCounts"]

sashimi = expand("rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
    rank=rank,
    as_type=as_type,
    jc_type=jc_type)

sashimi_short = expand("rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_index/{as_type}.event.list.txt",
    rank = rank,
    as_type=as_type,
    jc_type=jc_type)

bed12 = expand("rmats_out/MATS_output/{as_type}.MATS.{jc_type}.bed",
        as_type=as_type,
        jc_type=jc_type)

rule all:
    input: sashimi_short + bed12 

ruleorder: rmats > parse > sashimiplot 

rule rmats:
    """
    Run rMATS to get Alternative splicing information
    """
    input:
        gtf = "/data/athersh/Danio_rerio.GRCz11.94.gtf",
        b1 = get_treatment_bams(),
        b2 = get_control_bams(),
    output:
        "rmats_out/MATS_output/{as_type}.MATS.{jc_type}.txt"
    params:
        length = 101,
        reading = "single"
    wrapper:
        "file:wrappers/rmats"

rule parse:
    """
    Get the AS events based on rank
    """
    input:
        "rmats_out/MATS_output/{as_type}.MATS.{jc_type}.txt",
    output:
        pvalue_bottom5 = "rmats_out/MATS_output/pvalue_bottom5/{as_type}.MATS.{jc_type}.txt",
        fdr_bottom5 = "rmats_out/MATS_output/fdr_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_bottom5 = "rmats_out/MATS_output/incleveldifference_bottom5/{as_type}.MATS.{jc_type}.txt",
        incleveldifference_top5 = "rmats_out/MATS_output/incleveldifference_top5/{as_type}.MATS.{jc_type}.txt",
    wrapper:
        "file:wrappers/parse"

rule sashimiplot:
    """
    Create the Sashimiplots to visualize the AS events
    """
    input: 
        rmats = "rmats_out/MATS_output/{rank}/{as_type}.MATS.{jc_type}.txt",
        b1 = get_treatment_bams(),
        b2 = get_control_bams()
    output:
        dir = directory("rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot"),
#        subdir = directory("rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/"),
        event = "rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_index/{as_type}.event.list.txt"
    params:
        as_ =  "{as_type}"
#    log:
#        "logs/sashimiplot/{sample}/{rank}/{as_type}/{jc_type}"
    wrapper:
        "file:wrappers/sashimiplot"

#rule pdfunite:
#    """
#    Combine the pdf's into a single file
#    """
#    input:
#        directory("rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/"),
#    output: 
#        "rmats_out/Sashimi_plots/{rank}/{as_type}.MATS.{jc_type}/Sashimi_plot/combined.pdf",
#    log:
#        "logs/pdfunite/{rank}/{as_type}/{jc_type}"
#    shell:
#        "pdfunite {input}* {output} &> {log}"

rule convert_to_bed12:
    """
    For each output combination from rMATS, convert the events file to bed12 format
    """
    input: "rmats_out/MATS_output/{as_type}.MATS.{jc_type}.txt"
    output: "rmats_out/MATS_output/{as_type}.MATS.{jc_type}.bed"
    wrapper: "file:wrappers/convert_to_bed12"


# vim: ft=python
