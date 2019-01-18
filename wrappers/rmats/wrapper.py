from snakemake import shell

treatment = ",".join(snakemake.input.b1)

control = ",".join(snakemake.input.b2)

outdir = "/".join(str(snakemake.output).split("/")[:2])

shell("python2.7 rMATS.3.2.5/RNASeq-MATS.py "
    "-b1 {treatment} -b2 {control} "
    "-t {snakemake.params.reading} "
    "-len {snakemake.params.length} "
    "-gtf {snakemake.input.gtf} "
    "-o {outdir} &> {snakemake.log}")
