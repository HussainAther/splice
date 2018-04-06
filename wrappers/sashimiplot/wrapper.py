from snakemake import shell

treatment = ",".join(snakemake.input.b1)

control = ",".join(snakemake.input.b2)

shell("python rmats2sashimiplot/src/rmats2sashimiplot/rmats2sashimiplot.py "
      "--b1 {treatment} --b2 {control} "
      "--l1 \"experiment\" --l2 \"control\" "
      "-o {snakemake.output.dir} -e {snakemake.input.rmats} "
      "-t {snakemake.params.as_} &> {snakemake.log}")

shell("touch {snakemake.output.event}")
