#!/usr/bin/python
import os
import shutil
import datetime
import pandas as pd
import numpy as np
from pathlib import Path
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator.scripts import CNVMatrixGenerator as scna
from collections import Counter
import random


FEATURES_LIST = ['NCTG', 'NCGT', 'DEL_5_MH', 'LOH.1.40Mb', '3-9:HET.10.40Mb', '2-4:HET.40Mb']

def features(genome='GRCh38', exome=True, SNV_DIR=None,INDELS_DIR=None,CNV_DIR=None,RESULT_DIR=None,cnv_file_type ='ASCAT', bootstrap=True,nreplicates=20):

    '''
        Extracts all features to determine HRD status for samples
        Required input parameters:
            1. genome: Genome build of VCF files. Currently, the pipeline only supports GRCh37 and GRCh38 (default: GRCh37)
            2. exome: Is the input data exome or not (default: True)
            3. SNV_DIR: Directory path to SNV files. If the files contain both indels and SNVs, then set INDEL_DIR=None.
            4. INDELS_DIR: Directory path to indel files. If SNV_DIR contains mutation files that contain both indels and SNVs, then do not provide a path to this directory. If the indel calls and SNV calls are separate, then provide the path to this directory.
            5. CNV_DIR: Directory path to copy number files.
            6. RESULT_DIR: Path to the output directory where the output and the log files for HRProflier will be available.
            7. cnv_file_type: File type for CNV files provided. Currently, the pipeline supports types ASCAT, ASCAT-NGS, SEQUENZA, and ABSOLUTE.
        Optional parameters:
            8. bootstrap: Simulate features per sample based on the sample-weighted probability (default: False)
            9. nreplicates: The number of replicates per sample (default: 20)

        Output:
            Pandas dataframe with the following feature columns:
                1. N[C>T]G: Proportion of C:G>T:A single base substitutions at 5’-NpCpG-3’ context (N[C>T]G)
                2. N[C>G]T: Proportion of C:G>G:C single base substitutions at 5’-NpCpT-3’ context (N[C>G]T)
                3. DEL:5:MH: Proportion of deletions spanning at least 5bp at microhomologies (abbreviated as DEL.5.MH). For exome samples, the number of deletions spanning at least 5bp at microhomologies are extracted
                4. LOH:1-40Mb: Proportion of genomic segments with loss of heterozygosity (LOH) with sizes between 1 and 40 megabases (LOH:1-40Mb)
                5. 2-4:HET:>40Mb: Proportion of heterozygous genomic segments with TCN between 2 and 4 and sizes above 40 megabases (2-4:Het:>40Mb)
                6. 3-9:HET:10-40Mb: Proportion of heterozygous genomic segments with TCN between 3 and 9 and sizes between 10 and 40 megabases (3-9:HET:10-40Mb)

    '''

    print_progress('****************** Extracting ' +  "\x1B[3m" + 'HRProfiler' + "\x1B[3m" + ' Features ******************')


    if not check_input(SNV_DIR, INDELS_DIR, CNV_DIR):
        return None
        
    print_progress('All input files present!')
 
    sbs_df = extract_sbs_features(SNV_DIR,RESULT_DIR,exome=exome,genome=genome,bootstrap=bootstrap,nreplicates=nreplicates,check_for_indels=(INDELS_DIR is None  or  len(os.listdir(INDELS_DIR) ) == 0))
    indels_df = extract_indel_features(INDELS_DIR,RESULT_DIR,exome=exome,genome=genome,bootstrap=bootstrap,nreplicates=nreplicates) if INDELS_DIR else None
    cnv_df = extract_cnv_features(CNV_DIR, RESULT_DIR,cnv_file_type=cnv_file_type,bootstrap=bootstrap,nreplicates=nreplicates) if CNV_DIR else None

    if ( any(df is None for df in [sbs_df, cnv_df, indels_df]) and (INDELS_DIR is not None)) or (any(df is None for df in [sbs_df, cnv_df]) and (INDELS_DIR is None)):
        print_progress('Process terminiated! Either the SNV or the CNV features are not provided.')
        return None

    if indels_df is not None: 
        features_df = merge_features(sbs_df, cnv_df, indels_df, bootstrap) 
    else: 
        features_df = merge_features(sbs_df, cnv_df, None, bootstrap)
    print_progress('************************* Done extracting HRProfiler features *************************')
    return features_df



def check_input(SNV_DIR, INDELS_DIR, CNV_DIR):

    if (INDELS_DIR is None) and (SNV_DIR is not None  and  len(os.listdir(SNV_DIR) ) != 0):
        indels=None
        print_progress('Path to only SNVs directory is provided. Assuming the files in the directory contain both indels and SNVs.')
                
    elif (SNV_DIR is None or  len(os.listdir(SNV_DIR) ) == 0) and  (INDELS_DIR is None  or  len(os.listdir(INDELS_DIR) ) == 0):
        print_progress('Either input directory for both indels and SNV  is empty or the path provided is wrong.Please provide the path to SNV directory if VCF files  contain both indels and SNVs. Otherwise, path to both SNV and INDEL directories is required.')
        print_progress('Process terminated')
        return False
        
    if CNV_DIR is None  or  len(os.listdir(CNV_DIR) ) == 0:
        print_progress('Either input directory for copy number is empty or the path provided is wrong.')
        print_progress('Process terminated')
        return False
    return True


def merge_features(sbs_df, cnv_df, indel_df, bootstrap):

    if bootstrap:
        common_cols = ['samples', 'samples_iter']
    else:
        common_cols = ['samples']

    if sbs_df is None:
        return None
    elif indel_df is None:
        merged_df = pd.merge(sbs_df, cnv_df, on=common_cols)
    else:
        merged_df = pd.merge(sbs_df, pd.merge(indel_df, cnv_df, on=common_cols), on=common_cols)
    
    if bootstrap:
        return merged_df[['samples', 'samples_iter'] + FEATURES_LIST]
    else:
        return merged_df[['samples']  + FEATURES_LIST]


def extract_sbs_features(SNV_DIR,RESULT_DIR=None,exome=True,genome="GRCh38",bootstrap=False,nreplicates=20, bed_file=None, check_for_indels=False):


    """
    Extracts SNV HRProfiler Features.

    Parameters:
    - SNV_DIR (str): The directory containing the SNV data.
    - RESULT_DIR (str, optional): The output directory. Defaults to None.
    - exome (bool, optional): Specifies whether the data is exomic. Defaults to True.
    - genome (str, optional): Genome version. Defaults to "GRCh38".
    - bootstrap (bool, optional): Specifies whether to perform bootstrap. Defaults to False.
    - nreplicates (int, optional): Number of replicates for bootstrap. Defaults to 20.
    - bed_file (str, optional): The path to the exome BED file. Defaults to None.
    - check_for_indels (bool, optional): Specifies whether to check for indels. Defaults to False.

    Returns:
    - pandas.DataFrame: DataFrame containing extracted features.
    """
  

    print_progress('Extracting SNV HRProfiler Features:')

    try:
        matrices = matGen.SigProfilerMatrixGeneratorFunc('SNVs', genome, SNV_DIR, plot=False, exome=exome,bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False)
    except Exception as e:
        print_progress('SigProfilerMatrixGeneratorFunc Error.')
        print_progress(f'Error Details: {str(e)}')
        print_progress('Process terminated')
        cleanup(SNV_DIR,RESULT_DIR)
        return None

    if '96' not in matrices:
        print_progress( 'snvs not found.Process terminated')
        cleanup(SNV_DIR,RESULT_DIR)
        return None

    sbs_df = preprocess_matrix(os.path.join(SNV_DIR, 'output/SBS/SNVs.SBS96.exome' if exome else 'output/SBS/SNVs.SBS96.all'))
    sbs_df = extract_features(data=sbs_df, regex=">",features=['NCTG', 'NCGT'], exome=exome)
    print_progress('Successfully extracted SBS features for the following samples: ' + " ".join(sbs_df['samples'].tolist()))


    if bootstrap:
        bootstrap_sbs_df = bootstrap_features(data=sbs_df, sample_col='samples', regex=">",nreplicates=nreplicates, features=['NCTG', 'NCGT'], exome=exome)
    else:
        sbs_df =  sbs_df.loc[:,['samples', 'NCTG', 'NCGT']]
    

    if 'ID' in matrices and check_for_indels:
        print_progress( 'Indels present in the VCF files')
        indel_df = preprocess_matrix(os.path.join(SNV_DIR, 'output/ID/SNVs.ID83.exome' if exome else 'output/ID/SNVs.ID83.all'))
        indel_df = extract_features(data=indel_df, regex='Del:|Ins:',features=['DEL_5_MH'], exome=exome)
        print_progress('Successfully extracted Indel feature for the following samples: ' + " ".join(indel_df['samples'].tolist()))

        if bootstrap:
            bootstrap_indel_df =  bootstrap_features(data=indel_df, sample_col='samples', regex="Del:|Ins:",nreplicates=nreplicates, features=['DEL_5_MH'], exome=exome)
        else:
            indel_df =  indel_df.loc[:,['samples','DEL_5_MH']]

        # add the indel features to the SBS dataframe:
        sbs_df = pd.merge(sbs_df, indel_df, on='samples', how='outer') 

    cleanup(SNV_DIR,RESULT_DIR)
    return pd.merge(bootstrap_sbs_df, bootstrap_indel_df, on='samples_iter') if ('96' in matrices) and ('ID' in matrices) and (check_for_indels) and bootstrap else (bootstrap_sbs_df if bootstrap else sbs_df)



def extract_indel_features(INDELS_DIR,RESULT_DIR=None,exome=True,genome="GRCh38",bootstrap=False,nreplicates=20):

    """
    Extracts Indel HRProfiler Features.

    Parameters:
    - INDELS_DIR (str): The directory containing the Indels data.
    - RESULT_DIR (str): The output directory. Defaults to None.
    - exome (bool,optional): Specifies whether the data is exomic. Defaults to False.
    - genome (str, optional): Genome version. Defaults to "GRCh38".
    - bootstrap (bool, optional): Specifies whether to perform bootstrap. Defaults to False.
    - nreplicates (int, optional): Number of replicates for bootstrap. Defaults to 20.

    Returns:
    - pandas.DataFrame: DataFrame containing extracted DEL_5_MH feature.
    """

    print_progress('Extracting Indel HRProfiler Features:')

    try:
        matrices = matGen.SigProfilerMatrixGeneratorFunc('Indels', genome, INDELS_DIR, plot=False, exome=exome,bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False,cushion=100)
    except Exception as e:
        print_progress('SigProfilerMatrixGeneratorFunc Error.')
        print_progress(f'Error Details: {str(e)}')
        print_progress('Process terminated')
        cleanup(INDELS_DIR,RESULT_DIR)
        return None

    if 'ID' not in matrices:
        print_progress( 'indels not found.Process terminated')
        cleanup(INDELS_DIR,RESULT_DIR)
        return None

    else:
        df = preprocess_matrix(os.path.join(INDELS_DIR, f"output/ID/Indels.ID83.{'exome' if exome else 'all'}"))
        df = extract_features(data=df, regex='Del:|Ins:',features=['DEL_5_MH'], exome=exome)
        print_progress('Successfully extracted Indel feature for the following samples: ' + " ".join(df['samples'].tolist()))

        #cleanup:
        cleanup(INDELS_DIR,RESULT_DIR)

        return bootstrap_features(data=df, sample_col='samples', regex="Del:|Ins:", nreplicates=nreplicates, features=['DEL_5_MH'], exome=exome) if bootstrap else df.loc[:, ['samples', 'DEL_5_MH']]




def extract_cnv_features(CNV_DIR=None, RESULT_DIR=None,cnv_file_type='ASCAT',bootstrap=False,nreplicates=20):

    """
    Extracts CNV HRProfiler Features.

    Parameters:
    - CNV_DIR (str, optional): The directory containing the CNV data. Defaults to None.
    - RESULT_DIR (str, optional): The output directory. Defaults to None.
    - cnv_file_type (str, optional): Type of CNV file. Defaults to 'ASCAT'.
    - bootstrap (bool, optional): Specifies whether to perform bootstrap. Defaults to False.
    - nreplicates (int, optional): Number of replicates for bootstrap. Defaults to 20.

    Returns:
    - pandas.DataFrame: DataFrame containing extracted CNV features.
    """

    print_progress('Extracting CNV HRProfiler Features:')
    os.system('awk ' + " 'FNR==1 && NR!=1{next;}{print}' " + CNV_DIR +'/*.txt' + ' > ' +  RESULT_DIR + '/combined.seg.txt')
    #os.system('sed -e ' + " '1s/Sample/sample/g' " +  RESULT_DIR + '/combined.seg.txt')
    df = pd.read_csv(RESULT_DIR + '/combined.seg.txt', sep="\t")
    df.rename(columns={ df.columns[0]: "sample" }, inplace = True)
    df.to_csv(RESULT_DIR + '/combined.seg.txt', sep="\t")

    try:
        scna.generateCNVMatrix(cnv_file_type, RESULT_DIR+'/combined.seg.txt', 'cnv', RESULT_DIR)
    except Exception as e:
        print_progress('CNVMatrixGenerator Error.Process terminated')
        print_progress(f'Error Details: {str(e)}')
        cleanup(INDELS_DIR,RESULT_DIR)
        return

    if os.path.exists(RESULT_DIR+'/cnv.CNV48.matrix.tsv'):
        df = preprocess_matrix(RESULT_DIR+'/cnv.CNV48.matrix.tsv')
        df = extract_features(data=df, regex='het:|LOH:|homdel:',features=['LOH.1.40Mb','3-9:HET.10.40Mb','2-4:HET.40Mb'], exome=False)
        print_progress('Sucessfully extracted CNV feature counts: '  + " ".join(df['samples'].tolist()))
        
        # cleanup
        remove(path=RESULT_DIR + '/combined.seg.txt', kind='file')
        remove(path=RESULT_DIR + '/cnv.CNV48.matrix.tsv', kind='file')

        return bootstrap_features(data=df, sample_col='samples', regex='het:|LOH:|homdel:', nreplicates=nreplicates, features=['LOH.1.40Mb', '3-9:HET.10.40Mb', '2-4:HET.40Mb']) if bootstrap else df.loc[:, ['samples', 'LOH.1.40Mb', '3-9:HET.10.40Mb', '2-4:HET.40Mb']]

############################################ HELPER FUNCTIONS #############################################

def extract_features(data=None, regex="het:|LOH:|homdel:",features=None, exome=True):

    # keep counts for exome microhomology-mediated deletions feature
    if 'DEL_5_MH' in features and exome:
        data['DEL_5_MH'] = data.filter(regex='5:Del:M:').sum(axis=1)

    #convert all the channels to proportions:
    data['tot']=data.loc[:,data.filter(regex=regex).columns.tolist()].sum(axis=1)
    data.loc[:,data.filter(regex=regex).columns.tolist()] = data.loc[:,data.filter(regex=regex).columns.tolist()].div(data.tot.where(data.tot != 0, np.nan)  , axis=0)
    
    for feature in features:
        if feature == 'DEL_5_MH' and not exome:
            data['DEL_5_MH'] = data.filter(regex='5:Del:M:').sum(axis=1)
        elif feature == 'LOH.1.40Mb':
            data['LOH.1.40Mb'] = np.sum(data.filter(regex='LOH:1Mb|LOH:10Mb|LOH:>40Mb'), axis=1)
        elif feature == '3-9:HET.10.40Mb':
            data['3-9:HET.10.40Mb'] = data.filter(regex='3-4:het:10Mb-40Mb|5-8:het:10Mb-40Mb|9\\+:het:10Mb-40Mb').sum(axis=1)
        elif feature == '2-4:HET.40Mb':
            data['2-4:HET.40Mb'] = data.filter(regex='2:het:>40Mb|3-4:het:>40Mb').sum(axis=1)
        elif feature == 'NCTG':
            data['NCTG'] = data.filter(regex='C>T]G').sum(axis=1)
        elif feature == 'NCGT':
            data['NCGT'] = data.filter(regex='C\\[C>G]T|T\\[C>G]T|G\\[C>G]T|A\\[C>G]T').sum(axis=1)
    return data


def bootstrap_features(data=None, sample_col='samples', regex='het:|LOH:|homdel:',nreplicates=20, features=None, exome=False):
    
    bootstrap_df = pd.DataFrame()
    
    for i in range(len(data[sample_col])):
        bootstrap_tmp = boostrap_per_sample(mutation_types=data.filter(regex=regex).columns.tolist(),weights=data.filter(regex=regex).iloc[i,:].tolist(), totCounts=int(data.filter(regex='tot').iloc[i,:].values[0]), nreplicates=nreplicates)
        bootstrap_tmp = extract_features(data=bootstrap_tmp, regex=regex,features=features, exome=exome)
        bootstrap_tmp[sample_col] = [data[sample_col].iloc[i]] * nreplicates
        bootstrap_tmp[sample_col+'_iter'] = bootstrap_tmp[sample_col] + "_" + bootstrap_tmp['iter']
        bootstrap_df = pd.concat([bootstrap_df, bootstrap_tmp.loc[:, ['samples', sample_col+'_iter'] + features]], axis=0)
        new_row = data.iloc[[i]]
        new_row[sample_col+'_iter'] = data[sample_col].iloc[i] + "_org"
        bootstrap_df = pd.concat([new_row.loc[:, [sample_col,sample_col+'_iter']+features] ,bootstrap_df.loc[:, [sample_col,sample_col+'_iter']+features]]).reset_index(drop=True)
    return bootstrap_df


def boostrap_per_sample(mutation_types=None, weights=None, totCounts=None, nreplicates=None):
    
    numMT = len(mutation_types)
    tmp = pd.DataFrame();
    tmp['index'] = [ i for i in range(numMT)]
    for i in range(nreplicates):
        tmp = pd.merge( tmp,
                      pd.DataFrame.from_dict(Counter(random.choices([ i for i in range(numMT)],
                      weights=weights,
                      k=totCounts)),orient='index').reset_index(),
                     on='index',how = 'outer')
        
    df = tmp.T
    df = df.iloc[1:df.shape[0], :]
    df.columns = mutation_types
    df['iter'] = [str(i) for i in range(nreplicates)]
    df.fillna(0, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def preprocess_matrix(filepath=None):
    if os.path.exists(filepath):
        df = pd.read_csv(filepath,sep="\t")
        df = df.T
        df.columns = df.iloc[0,:]
        df = df.iloc[1:df.shape[0], :]
        df['samples'] = df.index.tolist()
        df.reset_index(drop=True, inplace=True)
        return df
    else:
        print_progress('Error: 96 matrix path does not exist: ' + filepath)
        return None

def print_progress(msg):
    now = datetime.datetime.now()
    now_text = now.strftime("%Y-%m-%d %H:%M")
    text = "[{0}] {1} ".format(now_text, msg)
    print(text)


def move_logs(source_dir, destination_dir):
    logs_dir = os.path.join(source_dir, "logs")
    if os.path.isdir(logs_dir):
        logs_files = os.listdir(logs_dir)
        for file in logs_files:
            shutil.move(os.path.join(logs_dir, file), os.path.join(destination_dir, "logs", file))
        remove(logs_dir)

def remove(path,kind = 'dir'):
    if kind  == 'file':
        if os.path.exists(path):
            os.remove(path)
    elif kind  == 'dir':
        if os.path.isdir(path):
            shutil.rmtree(path)  # remove dir and all contains
    else:
        raise ValueError("file {} is not a file or dir.".format(path))


def cleanup(SNV_DIR,RESULT_DIR):
    remove(path=SNV_DIR+"/output", kind='dir')
    remove(path=SNV_DIR+"/input", kind='dir')
    move_logs(SNV_DIR,RESULT_DIR)


















