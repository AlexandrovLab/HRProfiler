# HRProfiler
**H**omologous **R**ecombination Profiler (HRProfiler) is a classification tool that predicts HRD status for  both whole-genome and exome-sequenced breast and ovarian samples using HRD-specific mutation and copy number features.

[![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) 
[![Build Status](https://app.travis-ci.com/AlexandrovLab/HRProfiler.svg?token=tnMyG42yezzqz5hbqp9x&branch=main)](https://app.travis-ci.com/AlexandrovLab/HRProfiler)
## Installation

```bash
# clone the development branch of the *HRProfiler* tool
git clone -b main git@github.com:AlexandrovLab/HRProfiler.git

# install the dependencies 
cd HRProfiler
python3 setup.py install
```

## Run HRProfiler
From within a python session, you can run HRProfiler as follows:
```python
$ python3
>> from HRProfiler.scripts import HRProfiler as HR
>> HR.HRProfiler(data_matrix,genome,exome,INDELS_DIR,SNV_DIR,CNV_DIR,RESULT_DIR,cnv_file_type,bootstrap,nreplicates,normalize,hrd_prop_thresh,plot_predictions,organ)
```
The layout of the parameters are as follows:
|Parameters           |Info           |
|:-------------|:-------------| 
| data_matrix  | An optional pandas dataframe with the following required feature columns (in the same order):<ul><li>N[C>T]G: Proportion of C:G>T:A single base substitutions at 5’-NpCpG-3’ context</li><li>N[C>G]T: Proportion of C:G>G:C single base substitutions at 5’-NpCpT-3’ context</li><li>DEL:5:MH: Proportion (WGS) or counts (WES) of deletions spanning at least 5bp at microhomologies</li><li>LOH:1-40Mb: Proportion of genomic segments with loss of heterozygosity (LOH) with sizes at least 1 megabase</li><li>2-4:HET:>40Mb: Proportion of heterozygous genomic segments with TCN between 2 and 4 and sizes above 40 megabases</li><li>3-9:HET:10-40Mb: Proportion of heterozygous genomic segments with TCN between 3 and 9 and sizes between 10 and 40 megabases</li></ul>|
| genome  | Genome build of snv and indel input files. Options include: GRCh37 and GRCh38 (default: GRCh38) | 
| exome  | Is the input data exome or not (default: True) | 
| SNV_DIR  | Directory path to snv vcf/tab-delimited files (default: None) | 
| INDELS_DIR  | Directory path to indel vcf/tab-delimited files. If SNV_DIR contains mutation files with both indels and SNVs, then set INDEL_DIR=None (default: None) | 
| CNV_DIR  | Directory path to allele-specific segmentation files (default: None) | 
| RESULT_DIR  | Path to the directory where HRProfiler will save the output and the log files (default: None) | 
| cnv_file_type  | File type for CNV files provided. Options include: 'ASCAT' , 'ASCAT-NGS', 'ABSOLUTE','BATTENBERG', 'SEQUENZA' (default: 'ASCAT') | 
| bootstrap  | Simulate features per sample based on the sample-weighted probability (default: False)  | 
| nreplicates  | Number of replicates to simulate per sample (default: 20)  |  
| normalize  | Normalize each feature column by pre-defined mean and standard deviation (default: True)  | 
| hrd_prop_thresh | HRD Probability threshold to classify a sample as HRD (default: 0.5) | 
| plot_predictions | plot a histogram with the HRD probability values for all samples (default: True) | 
| organ | Organ type for prediction. Options include 'BREAST' and 'OVARIAN' (default: 'BREAST')


## Examples

#### 1. Extract HRD features and predict HRD status for samples with input vcf and copy number files:

To determine HRD status for samples with snvs and indel vcf files and allele-specific copy number calls, we can run the following example code where input files have been provided for 5 breast WGS samples. 


```python
# Note: cd to parent directory that contains the HRProfiler folder before executing the command
from HRProfiler.scripts import HRProfiler as HR
HR.HRProfiler(data_matrix=None,
              genome='GRCh37', 
              exome=False, 
              INDELS_DIR='./HRProfiler/example/input/indels',
              SNV_DIR='./HRProfiler/example/input/mutations/',
              CNV_DIR='./HRProfiler/example/input/copynumber/', 
              RESULT_DIR='./HRProfiler/example/output_example1/',
              cnv_file_type = 'ASCAT',
              bootstrap=False, 
              nreplicates=20,
              normalize=True, 
              hrd_prob_thresh=0.5,
              plot_predictions=True,
              organ='BREAST')
```

#### 2. Predict HRD status for samples with preprocessed HRD features: 

To determine the HRD status for samples with pre-defined HRD features: ['NCTG', 'NCGT', 'DEL_5_MH', 'LOH.1.40Mb', '3-9:HET.10.40Mb','2-4:HET.40Mb'], run the following command: 

```python
# Note: cd to parent directory that contains the HRProfiler folder before executing the command.
from HRProfiler.scripts import HRProfiler as HR
import pandas as pd

data_matrix = pd.read_csv('./HRProfiler/example/input/example_data_matrix.txt', sep="\t")

HR.HRProfiler(data_matrix=data_matrix,
              genome='GRCh37', 
              exome=False, 
              INDELS_DIR=None,
              SNV_DIR=None,
              CNV_DIR=None, 
              RESULT_DIR='./HRProfiler/example/output_for_example2/',
              cnv_file_type='ASCAT',
              bootstrap=False, 
              nreplicates=20,
              normalize=True, 
              hrd_prob_thresh=0.5,
              plot_predictions=True,
              organ='BREAST')
```
##  HRProfiler Output 
HRProfiler generates a histogram with the HRD probabilities per sample and a tab-delimited table with the following columns:<br />
  |Columns        |Info           |
  |:------------- |:-------------| 
  | samples      | sample names  |
  | NCTG      | Proportion of C:G>T:A single base substitutions at 5’-NpCpG-3’ context      | 
  | NCGT | Proportion of C:G>G:C single base substitutions at 5’-NpCpT-3’ context       | 
  | DEL_5_MH | Proportion or total counts of deletions spanning at least 5bp at microhomologies      | 
  | LOH.1.40Mb      | Proportion of genomic segments with loss of heterozygosity (LOH) with sizes between 1 and 40 megabases |
  | 3-9:HET.10.40Mb     | Proportion of heterozygous genomic segments with TCN between 2 and 4 and sizes above 40 megabases    |
  | 2-4:HET.40Mb | Proportion of heterozygous genomic segments with TCN between 2 and 4 and sizes above 40 megabases |
  |hrd.prob | HRD Probability | 
  |prediction | HRD status | 

  Additional columns provided if bootstrap=True:
  |Columns        |Info          |
  |:------------- |:-------------| 
  |mean.hrd.prob | Average HRD probability across all replicates |
  |LCI.95.hrd.prob | Lower confidence interval (2.5%) HRD Probability | 
  |UCI.95.hrd.prob | Upper confidence interval (97.5%) HRD Probability | 

## License
Academic Software License: © 2024 University of California, San Diego (“Institution”). Academic or nonprofit researchers are permitted to use this Software (as defined below) subject to Paragraphs 1-4:

1.	Institution hereby grants to you free of charge, so long as you are an academic or nonprofit researcher, a nonexclusive license under Institution’s copyright ownership interest in this software and any derivative works made by you thereof (collectively, the “Software”) to use, copy, and make derivative works of the Software solely for educational or academic research purposes, and to distribute such Software free of charge to other academic or nonprofit researchers for their educational or academic research purposes, in all cases subject to the terms of this Academic Software License. Except as granted herein, all rights are reserved by Institution, including the right to pursue patent protection of the Software.

2.	Any distribution of copies of this Software -- including any derivative works made by you thereof -- must include a copy (including the copyright notice above), and be made subject to the terms, of this Academic Software License; failure by you to adhere to the requirements in Paragraphs 1 and 2 will result in immediate termination of the license granted to you pursuant to this Academic Software License effective as of the date you first used the Software.

3.	IN NO EVENT WILL INSTITUTION BE LIABLE TO ANY ENTITY OR PERSON FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF INSTITUTION HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. INSTITUTION SPECIFICALLY DISCLAIMS ANY AND ALL WARRANTIES, EXPRESS AND IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE IS PROVIDED “AS IS.” INSTITUTION HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS OF THIS SOFTWARE.

4.	Any academic or scholarly publication arising from the use of this Software or any derivative works thereof will include the following acknowledgment:  The Software used in this research was created by Alexandrov Lab of University of California, San Diego. © 2022 University of California, San Diego.
