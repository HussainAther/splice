from snakemake import shell
import pandas as pd

df = pd.DataFrame.from_csv(str(snakemake.input), sep="\t")

tsv = pd.DataFrame.to_csv(df.sort_values(["PValue"]).head(), sep="\t")

with open(snakemake.output.pvalue_bottom5, "w") as file:
    file.write(tsv)

tsv = pd.DataFrame.to_csv(df.sort_values(["FDR"]).head(), sep="\t")

with open(snakemake.output.fdr_bottom5, "w") as file:
    file.write(tsv)

tsv = pd.DataFrame.to_csv(df.sort_values(["IncLevelDifference"]).head(), sep="\t")

with open(snakemake.output.incleveldifference_bottom5, "w") as file:
    file.write(tsv)

tsv = pd.DataFrame.to_csv(df.sort_values(["IncLevelDifference"]).tail(), sep="\t")

with open(snakemake.output.incleveldifference_top5, "w") as file:
    file.write(tsv)
