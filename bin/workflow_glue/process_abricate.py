#!/usr/bin/env python
"""Functions to help us parse files for RSV Abricate typing results."""

################################################################################################################
# Last modified: 2025-03-04, Shania Wu
################################################################################################################

import json

import pandas as pd

import pyranges as pr

from .util import get_named_logger, wf_parser  # noqa: ABS101

def get_total_coverage(abricate):
    # Create Pyranges object from abricate hits
    abricate_hits = pr.PyRanges(
                        chromosomes=abricate['PRODUCT'].iloc[0], 
                        starts=abricate['START'].values, 
                        ends=abricate['END'].values)
    
    # Merge overlapping intervals and calculate total covered bases
    total_covered_bases = abricate_hits.merge(strand=False).length
    print(f"Total coverage of hits: {total_covered_bases} bp")
    
    # Extract the total length of the reference genome from the first hit
    reference_length = int(abricate['COVERAGE'].iloc[0].split('/')[-1])  # Assumes COVERAGE is in "start-end/length" format
    print(f"Reference length: {reference_length} bp")

    # Calculate true %COVERAGE 
    percent_coverage = (total_covered_bases / reference_length) * 100
    print(f"% coverage: {percent_coverage}%")
    
    return(round(percent_coverage,2))
    
        
def parse_typing_file(typing_file):
    """Summarise abricate results."""
    abricate = pd.read_csv(typing_file, delimiter="\t", keep_default_na=False)
    
    ## Initialize variables in case missed in checks
    type_as_A = False
    type_as_B = False
    
    ## check if typing_file empty
    if abricate.empty:
        print("Abricate typing file empty. Unable to determine RSV type. Please check draft.consensus.fasta.")
        result = {
            'type': "None"
        } 
        return [result, abricate] # end here
    
    ## Get df with only RSV_TYPE_A and RSV_TYPE_B
    abricate_A = abricate[abricate['PRODUCT'] == "RSV_TYPE_A"]
    abricate_B = abricate[abricate['PRODUCT'] == "RSV_TYPE_B"]
    
    ## Combine hits for each type if more than one hit
    if len(abricate_A) > 1: # more than one hit 
        print("More than one hit for A detected: " +  str(len(abricate_A)) + "\nCombining hits...")
        
        total_coverage_A = get_total_coverage(abricate_A)
        
        dict_A = {'#FILE': abricate_A['#FILE'].iloc[0], 
                    'SEQUENCE': abricate_A['SEQUENCE'].iloc[0], 
                    'START': abricate_A['START'].iloc[0],  # START of first hit
                    'END': abricate_A['END'].iloc[-1],  # END of last hit
                    'STRAND': abricate_A['STRAND'].iloc[0], 
                    'GENE': abricate_A['GENE'].iloc[0], 
                    'COVERAGE': 'combined_hits', 
                    'COVERAGE_MAP': 'NA', 
                    'GAPS': 'NA', 
                    '%COVERAGE': total_coverage_A,  # true coverage calcuated in get_total_coverage
                    '%IDENTITY': round(abricate_A['%IDENTITY'].mean(),2),  # avg of all hits for type (biased)
                    'DATABASE': abricate_A['DATABASE'].iloc[0], 
                    'ACCESSION': abricate_A['ACCESSION'].iloc[0], 
                    'PRODUCT': abricate_A['PRODUCT'].iloc[0], 
                    'RESISTANCE': abricate_A['RESISTANCE'].iloc[0]}
        abricate_A_grouped = pd.DataFrame([dict_A])
    else: # one hit or no hits
        abricate_A_grouped = abricate_A
        
    if len(abricate_B) > 1: # more than one hit 
        print("More than one hit for B detected: " +  str(len(abricate_B)) + "\nCombining hits...")
        
        total_coverage_B = get_total_coverage(abricate_B)

        dict_B = {'#FILE': abricate_B['#FILE'].iloc[0], 
                    'SEQUENCE': abricate_B['SEQUENCE'].iloc[0], 
                    'START': abricate_B['START'].iloc[0],  # START of first hit
                    'END': abricate_B['END'].iloc[-1],  # END of last hit
                    'STRAND': abricate_B['STRAND'].iloc[0], 
                    'GENE': abricate_B['GENE'].iloc[0], 
                    'COVERAGE': 'combined_hits', 
                    'COVERAGE_MAP': 'NA', 
                    'GAPS': 'NA', 
                    '%COVERAGE': total_coverage_B,  # true coverage calcuated in get_total_coverage
                    '%IDENTITY': round(abricate_B['%IDENTITY'].mean(),2),  # avg of all hits for type (biased)
                    'DATABASE': abricate_B['DATABASE'].iloc[0], 
                    'ACCESSION': abricate_B['ACCESSION'].iloc[0], 
                    'PRODUCT': abricate_B['PRODUCT'].iloc[0], 
                    'RESISTANCE': 'NA'}
        abricate_B_grouped = pd.DataFrame([dict_B])
    else: # one hit or no hits
        abricate_B_grouped = abricate_B
        
    ## Concat grouped A and B results (regardless or coverage) into one df for TSV output
    abricate_grouped = pd.concat([abricate_A_grouped, abricate_B_grouped])


    #### Check if total % genome coverage > 50% for each genome
    if len(abricate_A_grouped) > 0:
        A_total_cov = abricate_A_grouped['%COVERAGE'].iloc[0]
        if A_total_cov > 50:
            print(f"Genome coverage for RSV_TYPE_A combined hits > 50% ({A_total_cov}%), A_above_50 = True")
            A_above_50 = True
        else:
            A_above_50 = False
            print(f"Genome coverage for RSV_TYPE_A combined hits <= 50% ({A_total_cov}%), A_above_50 = False")
    else:
        A_above_50 = False 
        print("No type A hits.")
    
            
    if len(abricate_B_grouped) > 0:
        B_total_cov = abricate_B_grouped['%COVERAGE'].iloc[0]
        if B_total_cov > 50:
            print(f"Genome coverage for RSV_TYPE_B combined hits > 50% ({B_total_cov}%), B_above_50 = True")
            B_above_50 = True
        else:
            B_above_50 = False
            print(f"Genome coverage for RSV_TYPE_B combined hits <= 50% ({B_total_cov}%), B_above_50 = False")
    else:
        B_above_50 = False 
        print("No type B hits.")
   
    #### Further assign types based on A_above_50 and B_above_50:
    ## 1) If both A and B have > 50% coverage, type as mixedAB
    if A_above_50 and B_above_50:
        print(f"Both type A and B detected with % genome coverage > 50%. Typing as 'mixedAB'.")                
        type_as_A = True
        type_as_B = True
        
    ## 2) If only A has % genome coverage > 50%, type as A, regardless of whether B has hits <= 50% coverage
    elif A_above_50 and not B_above_50:
        type_as_A = True
        type_as_B = False
        print(f"Only genome coverage for RSV_TYPE_A combined hits > 50% ({A_total_cov}%), typing as 'RSV_TYPE_A'.")
    
    ## 3) If only B has % genome coverage > 50%, type as B, regardless of whether A has hits <= 50% coverage
    elif B_above_50 and not A_above_50:
        type_as_A = False
        type_as_B = True
        print(f"Only genome coverage for RSV_TYPE_B combined hits > 50% ({B_total_cov}%), typing as 'RSV_TYPE_B'.")
        
    ## 4) If A and B both dont have > 50% coverage, do further checks
    elif not A_above_50 and not B_above_50:
        
        # there are hits <= 50% coverage for both A and B => mixedAB
        if len(abricate_A_grouped) > 0 and len(abricate_B_grouped) > 0: 
            print(f"Both RSV_TYPE_A ({A_total_cov}) and RSV_TYPE_B ({B_total_cov}) have hits but with <= 50% genome coverage, assigning as 'mixedAB'.")
            type_as_A = True
            type_as_B = True
    
        # only A has hits <= 50% coverage; B has no hits
        elif len(abricate_A_grouped) > 0 and len(abricate_B_grouped) == 0:
            print("Only RSV_TYPE_A has hits but with <= 50% genome coverage, typing as 'RSV_TYPE_A'.")
            type_as_A = True
            type_as_B = False
        
        # only B has hits <= 50% coverage; A has no hits
        elif len(abricate_A_grouped) == 0 and len(abricate_B_grouped) > 0:
            print("Only RSV_TYPE_B has hits but with <= 50% genome coverage, typing as 'RSV_TYPE_B'.")
            type_as_A = False
            type_as_B = True
        
        # no hits for A and B at any coverage (should have been caught by previous if abricate.empty check but just in case)
        else:
            print("No hits at any coverage for RSV_TYPE_A or RSV_TYPE_B, typing as 'None'.")
            type_as_A = False
            type_as_B = False
            

    ## Assign result for JSON output based on type_as_A and type_as_B
    if type_as_B and type_as_A: ## edge case if both types are detected and > 50% coverage
        print("Both type A and B detected with % genome coverage > 50%. Typing as 'mixedAB'.")                
        result = {
            'type': "mixedAB"
        }
    elif type_as_A and not type_as_B: # type A only
        print(f"Typing result: RSV_TYPE_A") 
        result = {
            'type': "RSV_TYPE_A"
        }
        
    elif type_as_B and not type_as_A: # type B only
        print(f"Typing result: RSV_TYPE_B") 
        result = {
            'type': "RSV_TYPE_B"
        }
    else:
        print("Unable to determine RSV type. Please check draft.consensus.fasta.")
        result = {
            'type': "None"
        } 
    return [result, abricate_grouped]


def main(args):
    """Run the entry point."""
    logger = get_named_logger("process_abricate")
    result_json = parse_typing_file(args.typing)[0]
    result_tsv = parse_typing_file(args.typing)[1]

    with open(args.output_json, 'w') as f:
        f.write(json.dumps(result_json, indent=4))

    result_tsv.to_csv(args.output_tsv, sep="\t", index=False)

    logger.info(f"Typing result written to {args.output_json} and {args.output_tsv}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("process_abricate")
    parser.add_argument(
        "--typing", default=None,
        help="abricate typing results.")
    parser.add_argument(
        "--output_json", default=None,
        help="processed abricate results JSON.")
    parser.add_argument(
        "--output_tsv", default=None,
        help="processed abricate results full TSV.")
    return parser