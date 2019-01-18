from snakemake import shell

treatment = ",".join(snakemake.input.b1)

control = ",".join(snakemake.input.b2)

shell("""rmats2sashimiplot 
         --b1 {treatment} --b2 {control} 
         --l1 \"experiment\" --l2 \"control\" 
         -o {snakemake.params.dir} -e {snakemake.input.rmats} 
         -t {snakemake.params.as_} &> {snakemake.log}""")
