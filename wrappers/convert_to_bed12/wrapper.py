import pandas as pd
import sys

input_file = str(snakemake.input)

df_o = pd.DataFrame()

as_type = snakemake.wildcards.as_type 
jc_type = snakemake.wildcards.jc_type

df = pd.read_table(input_file)
"""Sample rMATS output (SE)
ID	GeneID	geneSymbol	chr	strand	exonStart_0base	exonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
3568	"FBgn0052521"	NA	chrX	-	21223579	21223647	21211686	21211799	21279181	21279186	3568	2,1,2	4,3,3	1,3,2	0,0,0	100	50	7.512768185335972e-12	1.76099286264e-08	0.2,0.143,0.25	1.0,1.0,1.0	-0.802
660	"FBgn0031589"	NA	chr2L	-	4192632	4192968	4190558	4190824	4193309	4193481	660	3,3,4	0,0,0	1,3,1	3,3,2	100	50	1.854805198320264e-11	2.17383169243e-08	1.0,1.0,1.0	0.143,0.333,0.2	0.775
"""

"""(RI)
ID	GeneID	geneSymbol	chr	strand	riExonStart_0base	riExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IJC_SAMPLE_1	SJC_SAMPLE_1	IJC_SAMPLE_2	SJC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
889	"FBgn0015521"	NA	chr2L	+	2856376	2856705	2856376	2856568	2856630	2856705	889	2,2,2	5,3,1	5,2,2	0,0,0	100	50	5.945577363775101e-12	4.41161840392e-09	0.167,0.25,0.5	1.0,1.0,1.0	-0.694
"""

"""(A5SS)
ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IC_SAMPLE_1	SC_SAMPLE_1	IC_SAMPLE_2	SC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
2069	"FBgn0031453"	NA	chr2L	+	2752799	2752840	2752799	2752837	2753371	2753660	2069	12,10,10	11,12,6	0,0,0	19,19,19	53	50	1.2212453270876722e-15	1.60471635979e-12	0.507,0.44,0.611	0.0,0.0,0.0	0.519
2063	"FBgn0027835"	NA	chr2R	+	18412026	18412356	18412026	18412056	18413001	18413077	2063	7,8,4	3,2,2	9,11,6	0,0,0	350	50	1.647506575608304e-10	1.08241182017e-07	0.25,0.364,0.222	1.0,1.0,1.0	-0.721
"""

"""(A3SS)
ID	GeneID	geneSymbol	chr	strand	longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE	ID	IC_SAMPLE_1	SC_SAMPLE_1	IC_SAMPLE_2	SC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1056	"FBgn0038740"	NA	chr3R	-	19874344	19874527	19874344	19874510	19875397	19875678	1056	1,0,0	4,6,3	3,3,2	2,2,2	67	50	5.348690050199778e-06	0.00270694698382	0.157,0.0,0.0	0.528,0.528,0.427	-0.442
1487	"FBgn0085447"	NA	chr3L	+	5739443	5739726	5739464	5739726	5739105	5739334	1487	46,47,51	8,4,1	50,50,44	0,0,0	71	50	8.521554269558251e-06	0.00270694698382	0.802,0.892,0.973	1.0,1.0,1.0	-0.111
"""

"""(MXE)
ID	GeneID	geneSymbol	chr	strand	1stExonStart_0base	1stExonEnd	2ndExonStart_0base	2ndExonEnd	upstreamES	upstreamEE	downstreamES	downstreamEE	ID	IC_SAMPLE_1	SC_SAMPLE_1	IC_SAMPLE_2	SC_SAMPLE_2	IncFormLen	SkipFormLen	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
1551	"FBgn0033159"	NA	chr2R	-	7350580	7350704	7356106	7356230	7348130	7348975	7361717	7361867	1551	0,0,0	2,1,2	2,2,1	0,0,0	174	174	4.615197113366776e-13	1.05134190242e-09	0.0,0.0,0.0	1.0,1.0,1.0	-1.0
1593	"FBgn0033159"	NA	chr2R	-	7350157	7350281	7356106	7356230	7348130	7348975	7361717	7361867	1593	0,0,0	1,1,2	2,2,1	1,0,0	174	174	1.7155364906429327e-09	1.33867648537e-06	0.0,0.0,0.0	0.667,1.0,1.0	-0.889
"""
df_o["chrom"] = df["chr"] # chromosome
df_o["name"] = df["ID"].astype(str) + "." + as_type + "." + jc_type # unique identifier
df_o["score"] = 0  # score of 0
df_o["strand"] = df["strand"] # strand
df_o["itemRgb"] =  "255,0,0" 
df_o["FDR"] = df["FDR"]

as_type = snakemake.wildcards.as_type # Alternative splicing type

# For each alternative-splicing type, 

if as_type == "SE":
    df_o["chromStart"] = df["upstreamES"]
    df_o["chromEnd"] = df["downstreamEE"]
    df_o["thickStart"] = df["exonStart_0base"]
    df_o["thickEnd"] = df["exonEnd"]
    df_o["blockCount"] = 3 
    df_o["blockSizes"] = (df["upstreamEE"] - df["upstreamES"]).astype(str) + "," + (df["exonEnd"] - df["exonStart_0base"]).astype(str) + "," + (df["downstreamEE"] - df["downstreamES"]).astype(str) 
    df_o["blockStarts"] = (df["upstreamES"] - df_o["chromStart"]).astype(str) + "," + (df["exonStart_0base"] - df_o["chromStart"]).astype(str) + "," + (df["downstreamES"] - df_o["chromStart"]).astype(str)  

if as_type == "RI":
    df_o["chromStart"] = df["upstreamES"]
    df_o["chromEnd"] = df["downstreamEE"]
    df_o["thickStart"] = df["riExonStart_0base"]
    df_o["thickEnd"] = df["riExonEnd"] 
    df_o["blockCount"] = 3 
    df_o["blockSizes"] = (df["upstreamEE"] - df["upstreamES"]).astype(str) + "," + (df["downstreamEE"] - df["downstreamES"]).astype(str) + "," + (df["riExonEnd"] - df["riExonStart_0base"]).astype(str)
    df_o["blockStarts"] = (df["upstreamES"] - df_o["chromStart"]).astype(str) + "," + (df["downstreamES"] - df_o["chromStart"]).astype(str) + "," + (df["riExonStart_0base"] - df_o["chromStart"]).astype(str)

if (as_type == "A5SS" or as_type == "A3SS"):
    chromStart = []
    chromEnd = []
    thickStart = []
    thickEnd = []
    blockCount = []
    blockSizes = []
    blockStarts = []
    for index, row in df.iterrows():
        if row["shortES"] < row["flankingEE"]:
            chromStart.append(row["shortES"])
            chromEnd.append(row["flankingEE"])
            blockSizes.append(str(row["shortEE"] - row["shortES"]) +"," +str(row["longExonEnd"] - row["longExonStart_0base"]) + "," + str(row["flankingEE"] - row["flankingES"]))
            blockStarts.append(str(0) + "," + str(row["longExonStart_0base"] - chromStart[-1]) + "," + str(row["flankingES"] - chromStart[-1]))
            thickStart.append(row["longExonStart_0base"])
            thickEnd.append(row["shortEE"])
            blockCount.append(3)
        else:
            chromEnd.append(row["shortEE"])
            chromStart.append(row["flankingES"])
            blockSizes.append(str(row["flankingEE"] - row["flankingES"]) + "," + str(row["shortEE"] - row["shortES"]) + "," + str(row["longExonEnd"] - row["longExonStart_0base"]))
            blockStarts.append(str(0) + "," + str(row["shortES"] - chromStart[-1]) + "," + str(row["longExonStart_0base"] - chromStart[-1]))
            thickStart.append(row["longExonStart_0base"])
            thickEnd.append(row["shortEE"])
            blockCount.append(3)
    df_o["chromStart"] = chromStart
    df_o["chromEnd"] = chromEnd
    df_o["thickStart"] = thickStart
    df_o["thickEnd"] = thickEnd
    df_o["blockCount"] = blockCount
    df_o["blockSizes"] = blockSizes
    df_o["blockStarts"] = blockStarts

if as_type == "MXE":
    df_o["chromStart"] = df["upstreamES"]
    df_o["chromEnd"] = df["downstreamEE"]
    df_o["thickStart"] = df["1stExonStart_0base"]
    df_o["thickEnd"] = df["2ndExonEnd"]
    df_o["blockCount"] = 4
    df_o["blockSizes"] = (df["upstreamEE"] - df["upstreamES"]).astype(str) + "," + (df["1stExonEnd"] - df["1stExonStart_0base"]).astype(str) + "," + (df["2ndExonEnd"] - df["2ndExonStart_0base"]).astype(str) + "," + (df["downstreamEE"] - df["downstreamES"]).astype(str)
    df_o["blockStarts"] = (df["upstreamES"] - df_o["chromStart"]).astype(str) + "," + (df["1stExonStart_0base"] - df_o["chromStart"]).astype(str) + "," + (df["2ndExonStart_0base"] - df_o["chromStart"]).astype(str) + "," + (df["downstreamES"] - df_o["chromStart"]).astype(str) 

column_order = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts", "FDR"]

df = pd.DataFrame(df_o, columns=column_order)
df = df.loc[df["FDR"] < 0.05] # Keep only those with FDR less than .05
df = df.drop(columns=["FDR"])
tsv = pd.DataFrame.to_csv(df, sep="\t", index=False, header=False)

with open(str(snakemake.output), "w") as file:
    file.write(tsv)
