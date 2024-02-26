# Import statements
import sys
import platform
import time
import datetime
import os
import numpy as np
import pandas as pd
from HRProfiler.scripts.features import *
from HRProfiler.scripts.predict import *


def HRProfiler(data_matrix=None, genome='GRCh38',exome=True, SNV_DIR=None, INDELS_DIR=None,CNV_DIR=None,RESULT_DIR=None,cnv_file_type ='ASCAT',bootstrap=False, nreplicates=50,normalize=True, hrd_prob_thresh=0.4,plot_predictions=True,organ='BREAST'):
    
    '''
        HRProfiler Pipeline:
    
        Extracts all features to determine HRD status for samples
        Required input parameters:
            1. genome: Genome build of VCF files. Currently, the pipeline only supports GRCh37 and GRCh38 (default: GRCh37)
            2. exome: Is the input data exome or not (default: True)
            3. SNV_DIR: Directory path to SNV files. If the files contain both indels and SNVs, then set INDEL_DIR=None.
            4. INDELS_DIR: Directory path to indel files. If SNV_DIR contains mutation files that contain both indels and SNVs, then do not provide a path to this directory. If the indel calls and SNV calls are separate, then provide the path to this directory.
            5. CNV_DIR: Directory path to copy number files.
            6. RESULT_DIR: Path to the output directory where the output and the log files for HRProflier will be available.
            7. cnv_file_type: File type for CNV files provided. Currently, the pipeline supports types ASCAT, ASCAT-NGS, SEQUENZA, and ABSOLUTE.
    
        Optional Input parameters:
            1. data_matrix: Provide pandas data frame as input with the following feature columns:
                1. N[C>T]G: Proportion of C:G>T:A single base substitutions at 5’-NpCpG-3’ context (N[C>T]G)
                2. N[C>G]T: Proportion of C:G>G:C single base substitutions at 5’-NpCpT-3’ context (N[C>G]T)
                3. DEL:5:MH: Proportion of deletions spanning at least 5bp at microhomologies (abbreviated as DEL.5.MH). For exome samples, the number of deletions spanning at least 5bp at microhomologies are extracted
                4. LOH:1-40Mb: Proportion of genomic segments with loss of heterozygosity (LOH) with sizes between 1 and 40 megabases (LOH:1-40Mb)
                5. 2-4:HET:>40Mb: Proportion of heterozygous genomic segments with TCN between 2 and 4 and sizes above 40 megabases (2-4:Het:>40Mb)
                6. 3-9:HET:10-40Mb: Proportion of heterozygous genomic segments with TCN between 3 and 9 and sizes between 10 and 40 megabases (3-9:HET:10-40Mb)
    
            2. bootstrap: Simulate features per sample based on the sample-weighted probability (default: False)
            3. nreplicates: The number of replicates per sample (default: 20)
            4. normalize: normalize the proportions by the mean and standard deviation of each column (default: True)
            5. hrd_prob_thresh: HRD Probability threshold to classify  a sample as HRD (default: 0.3)
            6. plot_predictions: plot a histogram with the HRD probability values for all samples (default: True)
            7. organ: cancer type to predict
    
    '''
    
    print_progress('****************** Starting HRProfiler Pipeline ******************')
    
    # save input parameters to a log file
    init(data_matrix=data_matrix, genome=genome,exome=exome, SNV_DIR=SNV_DIR, INDELS_DIR=INDELS_DIR,CNV_DIR=CNV_DIR,RESULT_DIR=RESULT_DIR,cnv_file_type =cnv_file_type,bootstrap=bootstrap, nreplicates=50,normalize=True, hrd_prob_thresh=hrd_prob_thresh,plot_predictions=plot_predictions)
    
    if data_matrix is None:
        print_progress('Data matrix is empty. Now checking if the input mutation and copy number files are provided.')
        data_matrix = features(genome=genome, exome=exome, SNV_DIR=SNV_DIR,INDELS_DIR=INDELS_DIR,CNV_DIR=CNV_DIR,RESULT_DIR=RESULT_DIR,cnv_file_type =cnv_file_type, bootstrap=bootstrap,nreplicates=nreplicates)
        predict(data=data_matrix,exome=exome,hrd_prob_thresh=hrd_prob_thresh,RESULT_DIR=RESULT_DIR,normalize=normalize,plot_predictions=plot_predictions, bootstrap=bootstrap,organ=organ)
    
    # data matrix is not None and contains all the required HRD features
    else:
        predict(data=data_matrix,exome=exome,hrd_prob_thresh=hrd_prob_thresh,RESULT_DIR=RESULT_DIR,normalize=normalize,plot_predictions=plot_predictions, bootstrap=bootstrap,organ=organ)
    
    print_progress('************************* DONE *************************')
    return;



def init(data_matrix=None, genome='GRCh38', exome=True, SNV_DIR=None, INDELS_DIR=None, CNV_DIR=None, RESULT_DIR=None, cnv_file_type='ASCAT', bootstrap=False, nreplicates=50, normalize=True, hrd_prob_thresh=0.4, plot_predictions=True):
    
    '''
    Initialize HRProfiler
    '''
    time_stamp = datetime.date.today()
    output_log_path = os.path.join(RESULT_DIR, "logs")
    output_results_path = os.path.join(RESULT_DIR, "output")
    os.makedirs(output_log_path, exist_ok=True)
    os.makedirs(output_results_path, exist_ok=True)
    
    error_file = os.path.join(output_log_path, f'HRProfiler_{time_stamp}.err')
    log_file = os.path.join(output_log_path, f'HRProfiler_{time_stamp}.out')
    
    # Open log files
    tempErr = sys.stderr
    sys.stderr = open(error_file, 'w')
    log_out = open(log_file, 'w')
    
    log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
    log_out.write("-------System Info-------\n")
    log_out.write("Operating System Name: "+ platform.uname()[0]+"\n"+"Nodename: "+ platform.uname()[1]+"\n"+"Release: "+ platform.uname()[2]+"\n"+"Version: "+ platform.uname()[3]+"\n")
    log_out.write("\n-------Python and Package Versions------- \n")
    log_out.write("Python Version: "+str(platform.sys.version_info.major)+"."+str(platform.sys.version_info.minor)+"."+str(platform.sys.version_info.micro)+"\n")
    log_out.write("pandas version: "+pd.__version__+"\n")
    log_out.write("numpy version: "+np.__version__+"\n")
    log_out.write("\n-------Vital Parameters Used for the execution -------\n")
    log_out.write("data_matrix: {}\nGenome: {}\nexome: {}\nSNV_DIR: {}\nINDELS_DIR: {}\nCNV_DIR: {}\nRESULT_DIR: {}\ncnv_file_type: {}\nbootstrap: {}\nnreplicates: {}\nnormalize: {}\nhrd_prob_thresh: {}\nplot_predictions: {}\n".format( 'False' if data_matrix is None else 'True', genome, str(exome), str(SNV_DIR), str(INDELS_DIR),str(CNV_DIR),str(RESULT_DIR),str(cnv_file_type), str(bootstrap), str(nreplicates),str(normalize), str(hrd_prob_thresh),str(plot_predictions)))
    log_out.write("\n-------Date and Time Data------- \n")
    tic = datetime.datetime.now()
    log_out.write("Date and Clock time when the execution started: "+str(tic)+"\n\n\n")
    log_out.write("-------Runtime Checkpoints------- \n")
    log_out.close()
    return;

