import pandas as pd
def get_tair_coordinates(gff, arabidopsis_gene = FALSE, arabidopsis_exon = FALSE, \
                        arabidopsis_three_prime_UTR = FALSE, \
                        arabidopsis_five_prime_UTR = FALSE):
    """
    _summary_
    arabidopsis gff clean and get coordinates a datascience 
    approach which uses the dataframe and a faster access to the 
    coordinates and provides the gene, exon, three_prime_UTR and 
    five_prime_UTR nested lists with the ID and the nested coordinates
    as list. You can download the gff from the TAIR: TAIR10_GFF3_genes.gff

    Arguments:
        gff -- _a tair gff _
        gene -- _if you want gene coordinates_
        exon -- _if you want exon coordinates_
        three_prime_UTR -- _if you want three_prime_UTR_
        five_prime_UTR -- _if you want five_prime_UTR_
    """
   
    tair = pd.read_csv(gff, sep = "\t")
    renaming_tair = tair.rename(columns={"Chr1":"Chromosome", "chromosome": "gene_type", "1": \
                            "Start", "30427671": "End", "..1": "Strand", "ID=Chr1;Name=Chr1":"Gene_ID"})
    all_prediction_type_start = renaming_tair["Start"]
    all_prediction_type_end = renaming_tair["End"]
    gene_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "gene").dropna()
    gene_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "gene").dropna()["Start"].to_list()))
    gene_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "gene").dropna()["End"].to_list()))
    gene_type_gene_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]].where(renaming_tair["gene_type"] == "gene").dropna() 
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))
    gene_type_gene_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "gene").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: \
                                                               n.replace("ID=", "")))["Gene_ID"].to_list()
    exon_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "exon").dropna()
    exon_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] \
                                                                                    == "exon").dropna()["Start"].to_list()))
    exon_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "exon").dropna()["End"].to_list()))
    exon_type_exon_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))
    exon_type_exon_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "exon").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))["Gene_ID"].to_list()
    three_prime_UTR_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "three_prime_UTR").dropna()
    three_prime_UTR_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "three_prime_UTR").dropna()["Start"].to_list()))
    three_prime_UTR_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "three_prime_UTR").dropna()["End"].to_list()))
    three_prime_UTR_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                      where(renaming_tair["gene_type"] == "three_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))
    three_prime_UTR_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "three_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))["Gene_ID"].to_list()
    five_prime_UTR_type = renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == "five_prime_UTR").dropna()
    five_prime_UTR_type_start = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "five_prime_UTR").dropna()["Start"].to_list()))
    five_prime_UTR_type_end = list(map(int,renaming_tair[["gene_type", "Start", "End"]].where(renaming_tair["gene_type"] == \
                                                                                    "five_prime_UTR").dropna()["End"].to_list()))
    five_prime_UTR_ID = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "five_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))
    five_prime_UTR_ID_AGI = pd.DataFrame(renaming_tair[["gene_type", "Start", "End","Gene_ID"]]. \
                                     where(renaming_tair["gene_type"] == "five_prime_UTR").dropna() \
                                      ["Gene_ID"].apply(lambda n : n.split(";")[0]).apply(lambda n : n.split(".")[0]) \
                                          .apply(lambda n: n.replace("Parent=", "")).apply(lambda n: n.replace("ID=", "")))["Gene_ID"].to_list()
    arabidopsis_gene = []
    for i in range(len(gene_type_start)):
        arabidopsis_gene.append([[gene_type_gene_ID_AGI[i]],[gene_type_start[i], gene_type_end[i]]])
    arabidopsis_exon = []
    for i in range(len(exon_type_start)):
        arabidopsis_exon.append([[exon_type_exon_ID_AGI[i], exon_type_start[i], exon_type_end[i]]])
    arabidopsis_three_prime_UTR = []
    for i in range(len(three_prime_UTR_type_start)):
        arabidopsis_three_prime_UTR.append([three_prime_UTR_ID_AGI[i], three_prime_UTR_type_start[i], three_prime_UTR_type_end[i]])
    arabidopsis_five_prime_UTR = []
    for i in range(len(five_prime_UTR_type_start)):
        arabidopsis_five_prime_UTR.append([[five_prime_UTR_ID_AGI[i],five_prime_UTR_type_start[i], five_prime_UTR_type_end[i]]])  
    if gff and arabidopsis_gene:
        return arabidopsis_gene
    if gff and arabidopsis_exon:
        return arabidopsis_exon
    if gff and arabidopsis_three_prime_UTR:
        return arabidopsis_three_prime_UTR
    if gff and arabidopsis_five_prime_UTR:
        return arabidopsis_five_prime_UTR                          
