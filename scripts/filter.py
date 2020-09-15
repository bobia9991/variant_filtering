#!/usr/bin/env python

import pandas as pd
import argparse
import numpy as np
import csv

def fix_clinvar_gene(df):

  df['ClinVar_GENE'] = df['ClinVar_GENEINFO'].str.split(":",n=1,expand=True)[0]
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar_CLNHGVS'] = np.nan
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar_CLNSIG'] = np.nan
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar_CLNREVSTAT'] = np.nan
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar_CLNDN'] = np.nan
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar_GENEINFO'] = np.nan 
  df.loc[df['ClinVar_GENE'] != df['SYMBOL'],'ClinVar'] = np.nan 
  del df['ClinVar_GENE']

  return df

def read_genelist(file_name):
    
    gene_list = [line.rstrip('\n') for line in open(file_name)]
    
    return gene_list

def filter_genelist(df, file_name):
    
    gene_list = read_genelist(file_name)
    df_filter = df.loc[ df['SYMBOL'].isin(gene_list) ]
    
    return df_filter

def filter_stringent(df,chrm_string):
    
    no_assert = ['no_assertion_criteria_provided']
    paths = ['Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic','Pathogenic&_other']
    impacts = ['MODERATE', 'HIGH']
    up_down_cons = ['upstream_gene_variant', 'downstream_gene_variant']

    assertion_cond = (~df['ClinVar_CLNREVSTAT'].isin(no_assert))    
    path_cond = (df['ClinVar_CLNSIG'].isin(paths))    
    af_cond = ((df['MAX_AF'].isnull()) | (df['MAX_AF'] < 0.01))    
    impact_cond = (df['IMPACT'].isin(impacts))
    chrm_cond =  (df['CHROM'] == chrm_string)
    not_chrm_cond = (~(df['CHROM'] == chrm_string))
    chrm_af_cond = (df['AlleleFreqH'] < 0.01)

    up_down_cond = (~df['Consequence'].isin(up_down_cons))

    df_filter = df.loc[ (not_chrm_cond & path_cond & assertion_cond & af_cond & impact_cond) | (chrm_cond & path_cond & assertion_cond & chrm_af_cond & impact_cond & up_down_cond) ]
    
    return df_filter

def filter_relaxed(df, chrm_string):
    
    no_assert = ['no_assertion_criteria_provided']
    paths = ['Pathogenic','Likely_pathogenic','Pathogenic/Likely_pathogenic','Pathogenic&_other']
    benigns = ['Benign','Likely_benign','Benign/Likely_benign','protective']
    impacts = ['MODERATE', 'HIGH']
    up_down_cons = ['upstream_gene_variant', 'downstream_gene_variant']

    assertion_cond = (~df['ClinVar_CLNREVSTAT'].isin(no_assert))
    path_cond = (df['ClinVar_CLNSIG'].isin(paths))
    af_cond = ((df['MAX_AF'].isnull()) | (df['MAX_AF'] < 0.01))
    impact_cond = (df['IMPACT'].isin(impacts))
    benign_cond = ((~df['ClinVar_CLNSIG'].isin(benigns)) | (df['ClinVar_CLNSIG'].isnull()))
    chrm_cond =  (df['CHROM'] == chrm_string)
    not_chrm_cond = (~(df['CHROM'] == chrm_string))
    chrm_af_cond = ((df['AlleleFreqH'] < 0.01) | (df['AlleleFreqH'].isnull()))
    up_down_cond = (~df['Consequence'].isin(up_down_cons))

    df_filter = df.loc[ ( path_cond & assertion_cond & up_down_cond ) | ( not_chrm_cond & af_cond & impact_cond & benign_cond ) | ( chrm_cond & up_down_cond & chrm_af_cond & impact_cond & benign_cond )]
    
    return df_filter

def remove_cols(df):
    
    #for gencode VEP annotation and hmtnote annotation
    cols_to_remove = ["Gene","Feature_type","Feature","BIOTYPE","EXON","INTRON","HGVSc","HGVSp",
                      "CDS_position","Codons","Existing_variation","DISTANCE",
                      "STRAND","FLAGS","SYMBOL_SOURCE","HGNC_ID","CANONICAL",
                      "SIFT","PolyPhen","SOURCE","AFR_AF","AMR_AF","EAS_AF","EUR_AF",
                      "SAS_AF","AA_AF","EA_AF","gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF",
                      "gnomAD_EAS_AF","gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF","MAX_AF_POPS",
                      "CLIN_SIG","SOMATIC","PHENO","CADD_PHRED","CADD_RAW"]
    
    for col in cols_to_remove:
        del df[col]
    
    df = df.rename(columns={"Allele": "VEP_Allele",
                            "AlleleFreqH" : 'Mt_AlleleFreqH',
                            "Consequence": "VEP_Consequence",
                            "IMPACT":"VEP_IMPACT",
                            "SYMBOL":"VEP_SYMBOL",
                            "Protein_position":"VEP_Protein_position",
                            "Amino_acids":"VEP_Amino_acids",
                            "MAX_AF":"VEP_MAX_AF",
                            "cDNA_position":"VEP_cDNA_position"})
    return df

def npat(df):
    pats = 0
    for col in df:
        if col.endswith("[AD]:GQ:GT"):
            pats += 1
                
    return pats

def npat_variant(df):
    
    pats_with_path = 0
    for col in df:
        if col.endswith("[AD]:GQ:GT"):
            patient = col.strip("_[AD]:GQ:GT")
            l = list(df.loc[~df[col].isnull()][col])
            if len(l) > 0:
                pats_with_path += 1
                
    return pats_with_path

def format_biobank(tsv_in):
    f = tsv_in
    
    df = pd.read_csv(f, sep="\t", dtype={'ClinVar':object, "MAX_AF": "float64", "AlleleFreqH": "float64"})
    df = df.replace(np.nan, '', regex=True)
    
    patients = []
    
    for col in df:
        if col.endswith("[AD]:GQ:GT"):
            patients.append(col)
    
    ddt = df.T.to_dict()
    
    first_fields = ['CHROM','POS','REF', 'ALT','FILTER']
    
    second_fields = ['VEP_Allele','VEP_SYMBOL', 'VEP_IMPACT','VEP_Consequence','VEP_Protein_position','VEP_MAX_AF',
              'VEP_Amino_acids', 'VEP_cDNA_position',
              'ClinVar','ClinVar_CLNSIG','ClinVar_CLNREVSTAT','ClinVar_CLNDN']
    
    header = [['PATIENT'],first_fields,['REF_COVERAGE','ALT_COVERAGE','GENOTYPE'],second_fields]
    header = [item for sublist in header for item in sublist]
            
    data_final = {}
    for patient in patients:
        data_final[patient] = []
        
    for variant in ddt:
        for patient in patients:
            split_str = str(ddt[variant][patient]).split(":")
            nfields = len(split_str)
            
            if (nfields > 1):
                cov = split_str[0].split(",")
                ref_cov = cov[0]
                alt_cov = cov[1]
                genotype = split_str[2]
                
                patient_number = patient.strip("_[AD]:GQ:GT")
                out_string = [patient_number]
                
                for field in first_fields:
                    out_string.append(str(ddt[variant][field]))
                    
                out_string.append(str(ref_cov))
                out_string.append(str(alt_cov))
                out_string.append(str(genotype))
                
                for field in second_fields:
                    out_string.append(str(ddt[variant][field]))
                
                data_final[patient].append(out_string)
                
    out_file = open(f.replace("tsv","biobank.tsv"),'w')
    tsv_out = csv.writer(out_file, delimiter='\t')
    tsv_out.writerow(header)
    
    for patient in data_final:
        for variant in data_final[patient]:
            tsv_out.writerow(variant)
    
    out_file.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-genelist", dest="genelist", help="input genelist", required=False)    
    parser.add_argument("-tsv", dest="tsv_in", help="input tsv file", required=True)
    parser.add_argument("-build", dest="build", help="genome build; b37 and GRCh38 are supported", required=True)
    args = parser.parse_args()

    if args.build == 'b37':
      chrm_string = 'MT'
    elif args.build == 'GRCh38':
      chrm_string = 'chrM'
    else:
      sys.exit("build not supported")

    if args.genelist is None:
      genelist = ""
    else:
      genelist = args.genelist

    return genelist,args.tsv_in,chrm_string

def run_filter(tsv_in, genelist, chrm_string):
     
    df = pd.read_csv(tsv_in ,sep="\t", dtype={'ClinVar':object,"MAX_AF": "float64"},low_memory=False)
    df['AlleleFreqH'] = df['AlleleFreqH'].replace(r'^\.', np.nan, regex=True)
    df = df.astype({'AlleleFreqH': "float64"})
  
    df_fix = fix_clinvar_gene(df)
    out_unfiltered = tsv_in.replace("tsv","unfiltered.tsv")
    df_fix.to_csv(out_unfiltered, sep="\t",index=False)

    if genelist != "":
      df = filter_genelist(df, genelist)
      out_tmp = tsv_in.replace("tsv","genelist.tsv")
      df.to_csv(out_tmp, sep="\t",index=False)
    else: 
      out_tmp=tsv_in

    out1 = out_tmp.replace("tsv","removecols.stringent-filter.tsv")
    out2 = out_tmp.replace("tsv","removecols.relaxed-filter.tsv")
    
    df1 = filter_stringent(df,chrm_string)
    df1 = remove_cols(df1) 
    df1.to_csv(out1, sep="\t",index=False)
    
    df2 = filter_relaxed(df,chrm_string)
    df2 = remove_cols(df2)
    df2.to_csv(out2, sep="\t",index=False)

    format_biobank(out1)
    format_biobank(out2)

# main
genelist, tsv_in, chrm_string  = parse_args()
run_filter(tsv_in, genelist, chrm_string)

