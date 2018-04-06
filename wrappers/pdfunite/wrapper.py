from snakemake import shell

shell("pdfunite {snakemake.input}* {snakemake.output}")
