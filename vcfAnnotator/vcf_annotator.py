#!/usr/bin/env python3
import pandas as pd
import numpy as np
import vcfAnnotator.functions as fcn
from tqdm import tqdm
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def main():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        description="Annotate variants from a vcf file and outputs csv")
    parser.add_argument("vcf", type=str, help='vcf input file')
    parser.add_argument("output_prefix", type=str, help='prefix for output csv')
    opts = parser.parse_args()
    input_file = opts.vcf
    out_file = opts.output_prefix

    ## take in the vcf file as a pandas dataframe
    vcfdat = pd.read_csv(input_file, sep="\t", skiprows=48)

    ## column containing seq coverage, no of reads supporting variant, % of reads supporting variant, etc
    ## 14th entry is depth of seq coverage (TC)
    ## 16th entry is number of reads supporting variant (TR)
    ## 1th entry is the percentage of reads supporting variant vs supporting ref reads (or min allele)

    infoCol = list(vcfdat["INFO"])
    ## get the number of variants
    numVars = len(vcfdat)
    ## array for depth of sequence coverage
    TCarray = np.zeros(numVars)
    ##array for major allele reads (reads supporting variant)
    TRarray = np.zeros(numVars)
    ## array for % of reads supporting variant (major allele frequency)
    majAFarray = np.zeros(numVars)
    ## array for minor allele frequency (%)
    minAFarray = np.zeros(numVars)
    ##lists for hgvs notation for major and minor allele
    hgvslist = ["n/a" for x in minAFarray]
    hgvslist_min = ["n/a" for x in minAFarray]
    ##lists to store type of variation
    type_of_var = ["n/a" for x in minAFarray]
    type_of_var_min = ["n/a" for x in minAFarray]
    ##lists to store gene symbols
    gene_symbols = ["n/a" for x in minAFarray]
    ##lists to store effects of variants
    effects = ["n/a" for x in minAFarray]
    effects_min = ["n/a" for x in minAFarray]
    ## terms necessary for hgvs conversion
    chromosome = np.array(vcfdat["#CHROM"])
    pos = np.array(vcfdat["POS"])
    ref = np.array(vcfdat["REF"], dtype=str)
    alt = np.array(vcfdat["ALT"], dtype=str)
    for i in tqdm(np.arange(0, numVars)):
        varInfo = infoCol[i].split(";")
        TC = varInfo[14] #coverage depth
        TR = varInfo[17] #no. of reads for variant
        VF = varInfo[1] #variant frequency
        ## collect TC
        TC = int(TC.split("=")[1])
        ##store TC
        TCarray[i] = TC
        ## collect no. of reads for major and minor allele
        TR = TR.split("=")[1]
        TR = TR.split(",")
        ## collect the variant frequencies
        VF = VF.split("=")[1]
        VF = VF.split(",")
        ## get chromosome position, etc for hgvs conversion
        chrom1 = chromosome[i]
        pos1 = pos[i]
        ref1 = ref[i]
        alt1 = alt[i]
        ## check for minor allele
        if len(TR) > 1:
            ## decide which is major variant
            if VF[0] >= VF[1]:
                TR1 = int(TR[0])
                # TR2 = int(TR[1])
                VF1 = float(VF[0]) * 100
                VF2 = float(VF[1]) * 100
                # split alt two sequence for each variant
                alt1, alt2 = alt1.split(",")
            else:
                TR1 = int(TR[1])
                # TR2 = int(TR[0])
                VF1 = float(VF[1]) * 100
                VF2 = float(VF[0]) * 100
                # split alt two sequence for each variant
                alt2, alt1 = alt1.split(",")
            # collect values
            TRarray[i] = TR1
            majAFarray[i] = VF1
            minAFarray[i] = VF2
            # get hgvs and type of variation for major allele
            hgvs, var = fcn.vcf2hgvs(chrom1, pos1, ref1, alt1)
            hgvslist[i] = hgvs
            type_of_var[i] = var
            # get the effect and gene symbol from Ensembl REST API
            effect, gene_symbol = fcn.getVarInfo(hgvs)
            gene_symbols[i] = gene_symbol
            effects[i] = effect
            # get hgvs and type of variation for minor allele
            hgvs, var = fcn.vcf2hgvs(chrom1, pos1, ref1, alt2)
            hgvslist_min[i] = hgvs
            type_of_var_min[i] = var
            # get the effect from Ensembl REST API
            effect, _ = fcn.getVarInfo(hgvs)
            effects_min[i] = effect
        else:
            # just store no of reads and variant frequency if there is no minor allele
            TRarray[i] = int(TR[0])
            majAFarray[i] = float(VF[0])*100
            # get hgvs and type of variation
            hgvs, var = fcn.vcf2hgvs(chrom1, pos1, ref1, alt1)
            hgvslist[i] = hgvs
            type_of_var[i] = var
            # get the effect and gene symbol from Ensembl REST API
            effect, gene_symbol = fcn.getVarInfo(hgvs)
            gene_symbols[i] = gene_symbol
            effects[i] = effect

    # make an output csv file
    outputDf = pd.DataFrame(list(zip(chromosome, pos, ref, alt, hgvslist,
                                     TCarray, TRarray, majAFarray, gene_symbols, type_of_var, effects,
                                     minAFarray, hgvslist_min, type_of_var_min, effects_min)),
                            columns = ['Chrom','Pos', 'Ref', 'Alt', 'hgvs of variant', 'Depth of coverage',
                                       'No. of reads for variant', '% reads supporting variant', 'gene', 'variation', 'effect',
                                       'min allele freq (%)', 'min allele hgvs', 'min allele variation', 'min allele effect'])
    outputDf.to_csv(out_file+"_annotated_vcf_data.csv", header=True, index=False)

if __name__ == "__main__":
    main()