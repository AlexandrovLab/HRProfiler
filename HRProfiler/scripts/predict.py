import os
import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from matplotlib.backends.backend_pdf import PdfPages
from HRProfiler.scripts.plotting import *

FEATURES_LIST = ['NCTG', 'NCGT', 'DEL_5_MH', 'LOH.1.40Mb', '3-9:HET.10.40Mb', '2-4:HET.40Mb']
ORGANS=['BREAST', 'OVARIAN']

def predict(data=None, sample_col='samples', exome=True, normalize=True, hrd_prob_thresh=0.5, bootstrap=False, organ='BREAST', plot_predictions=True, RESULT_DIR=None):

    '''
        Predict the HRD Probability for the given samples

        Required input parameters:
            1. a pandas data frame with the following feature columns
                1. samples:  sample names
                2. N[C>T]G: Proportion of C:G>T:A single base substitutions at 5’-NpCpG-3’ context (N[C>T]G)
                3. N[C>G]T: Proportion of C:G>G:C single base substitutions at 5’-NpCpT-3’ context (N[C>G]T)
                4. DEL:5:MH: Proportion of deletions spanning at least 5bp at microhomologies (abbreviated as DEL.5.MH). For exome samples, the number of deletions spanning at least 5bp at microhomologies are extracted
                5. LOH:1-40Mb: Proportion of genomic segments with loss of heterozygosity (LOH) with sizes atleast 1 megabase (LOH:1-40Mb)
                6. 2-4:HET:>40Mb: Proportion of heterozygous genomic segments with total copy number (TCN) between 2 and 4 and sizes above 40 megabases (2-4:Het:>40Mb)
                7. 3-9:HET:10-40Mb: Proportion of heterozygous genomic segments with TCN between 3 and 9 and sizes between 10 and 40 megabases (3-9:HET:10-40Mb)
            2. sample_col: name of the sample column with all the sample names (default: samples)
            3. exome: input data exome or not (default: True)
            4. normalize: scale the features by the mean and standard deviation
            5. hrd_prob_thresh: HRD Probability threshold to classify  a sample as HRD (default: 0.4)
            6. bootstrap: Simulate features per sample based on the sample-weighted probability (default: False).
            7. organ: Name of the tissue-specific model to be used (default: 'BREAST'). Only models for breast (organ='BREAST') and ovarian (organ='OVARIAN') are currently supported
            6. plot_predictions: plot a histogram with the HRD probability values for all samples (default: True)
            7. RESULT_DIR: full path of the directiory where to store the output from HRProfiler

    '''

    # Check if all input features are present
    if not all(feature in data.columns for feature in FEATURES_LIST):
        print(f"There is something wrong with the input features provided. Please check if all the columns are provided: {FEATURES_LIST}")
        return

    # Load saved trained model
    path = os.path.join(os.path.dirname(__file__), 'Models')
    model_type = "WES" if exome else "WGS"
    if organ not in ORGANS:
        print('Incorrect organ-type provided. Please choose either BREAST or OVARIAN for HRD prediction')
        return
    model_path = f"{path}/{model_type}/HRProfiler_{model_type}_{organ}.sav"
    print(f'Using {"exome" if exome else "whole-genome"} {organ} model for prediction...')
    model = pickle.load(open(model_path, 'rb'))

    # Select features
    X = data.loc[:, FEATURES_LIST]

    # Normalize features
    if normalize:
        if exome:
            if organ == 'BREAST':
                centers = [0.23260921, 0.05634558, 0.8125, 0.19864466, 0.09258245, 0.25361812]
                scales = [0.10628014, 0.04774947, 1.83969742, 0.11137158, 0.06988917, 0.1826042 ]

            elif organ == 'OVARIAN':
                centers = [0.16380864, 0.05846835, 2.01648352, 0.26525727, 0.12587575, 0.11636021]
                scales = [0.08496781, 0.03108982, 2.42123226, 0.09050332, 0.07369654, 0.07210832]

        else:
            if organ == 'BREAST':
                centers = [0.12248847, 0.06329762, 0.0942639, 0.15455065, 0.06581463, 0.07507382]
                scales=[0.05377641, 0.04025697, 0.11243681, 0.0957664, 0.06134753, 0.07734589]

    custom_scaler = CustomScaler(centers, scales)
    X = custom_scaler.transform(X)

    # Prediction probabilities
    probs = model.predict_proba(X)
    predictions = (probs[:, 1] >= hrd_prob_thresh).astype(int)
    df = pd.DataFrame(data={sample_col: data[sample_col],'hrd.prob': probs[:, 1],'prediction': predictions })
    df = pd.concat([df, data.loc[:, FEATURES_LIST]], axis=1)
    
    if bootstrap:
        df = process_bootstrap(df, data, sample_col, RESULT_DIR, hrd_prob_thresh,model_type,organ)

    df.to_csv(f'{RESULT_DIR}output/hrd.predictions.organ.{organ.lower()}.model_type.{model_type.lower()}.txt', sep='\t')
    print(f"Successfully predicted HRD probability for the following samples: {', '.join(df[sample_col].unique())}")

    if plot_predictions:
        plot_pred(df, 'mean.hrd.prob' if bootstrap else 'hrd.prob', sample_col, RESULT_DIR, bootstrap, hrd_prob_thresh, model_type,organ)
    return;

def process_bootstrap(df, data, sample_col, RESULT_DIR, hrd_prob_thresh,model_type,organ):
    pp = PdfPages(f'{RESULT_DIR}output/hrd.probability.organ.{organ.lower()}.model_type.{model_type.lower()}.bootstrap.pdf')
    df[sample_col + '_iter'] = data[sample_col +'_iter']

    updated_df = pd.DataFrame()
    for sample in df[sample_col].unique():
        tmp = df.loc[(df[sample_col] == sample) & ~(df[sample_col + '_iter'].str.contains(sample + '_org')), 'hrd.prob']
        org_prob = df.loc[(df[sample_col] == sample) & (df[sample_col + '_iter'].str.contains(sample + '_org')), 'hrd.prob'].values[0]

        # Confidence intervals
        alpha = 0.95
        LCI = max(0.0, np.percentile(tmp, (1.0 - alpha) / 2.0 * 100))
        UCI = min(1.0, np.percentile(tmp, (alpha + (1.0 - alpha) / 2.0) * 100))

        updated_df = pd.concat([updated_df, pd.DataFrame([[sample, org_prob, np.mean(tmp), LCI, UCI]])], axis=0)
        plot_hrd_prob_per_sample(pp=pp, data=tmp, title=sample, LCI=LCI, UCI=UCI, org_prob=org_prob)

    pp.close()
    updated_df.columns = [sample_col, 'hrd.prob', 'mean.hrd.prob', 'LCI.95.hrd.prob', 'UCI.95.hrd.prob']
    updated_df = pd.merge(data.loc[(data[sample_col + '_iter'].str.contains('_org')),:], updated_df, on=sample_col)
    updated_df['prediction'] = (updated_df['mean.hrd.prob'] >= hrd_prob_thresh).astype(int)
    return updated_df.loc[:, [sample_col] + ['hrd.prob', 'mean.hrd.prob', 'LCI.95.hrd.prob', 'UCI.95.hrd.prob', 'prediction'] + FEATURES_LIST]



class CustomScaler:
    def __init__(self, custom_means, custom_stds):
        self.mean_ = custom_means
        self.scale_ = custom_stds  # Use 'scale_' for standard deviation

    def transform(self, X):
        # Assuming X is a DataFrame
        return (X - self.mean_) / self.scale_
