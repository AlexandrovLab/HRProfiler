import os
from HRProfiler.scripts import HRProfiler as HR

# set up the input and the output paths
HRPROFILER_EXAMPLE_PATH = os.path.dirname(os.path.abspath(__file__))
input_path = os.path.join(HRPROFILER_EXAMPLE_PATH, "input/")
output_path = os.path.join(HRPROFILER_EXAMPLE_PATH, "output/")

# create the result directory if it does not exist
if not os.path.exists(output_path):
    os.makedirs(output_path)

HR.HRProfiler(data_matrix=None,
              genome='GRCh38', 
              exome=False, 
              INDELS_DIR=input_path + 'indels/',
              SNV_DIR=input_path + 'mutations/',
              CNV_DIR=input_path + 'copynumber/', 
              RESULT_DIR=output_path,
              cnv_file_type = 'ASCAT',
              bootstrap=False, 
              nreplicates=20,
              normalize=True, 
              hrd_prob_thresh=0.5,
              plot_predictions=True,
              organ='BREAST')
