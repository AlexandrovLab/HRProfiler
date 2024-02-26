# import statements:
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pylab
import scipy.stats as stats


def plot_pred(df=None, pred_col='HRD_PROB', samples_col='Samples',RESULT_DIR=None, bootstrap=False, hrd_prob_thresh=0.4,model_type='WES',organ='BREAST' ):
    fig=plt.figure()
    df.sort_values(by=[pred_col], ascending=False,inplace=True)
    plt.bar(df[samples_col].tolist(),df.loc[:,pred_col],color='gray')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.xlabel('Samples')
    plt.axhline(y = hrd_prob_thresh, color = 'k', linestyle = ':')
    if bootstrap:
        plt.ylabel('Average HRD Probability')
    else:
        plt.ylabel('HRD Probability')
    plt.xticks(size=6)
    plt.savefig(RESULT_DIR+'output/hrd.probability.organ.'+str(organ.lower())+'.model_type.'+str(model_type.lower())+'.pdf',bbox_inches='tight', dpi=300)


def plot_probQQ(pp=None, data=None ):
    fig=plt.figure(figsize=(2,2))
    stats.probplot(data, dist="norm", plot=pylab)
    pp.savefig(plt.gcf())


def plot_hrd_prob_per_sample(pp=None, data=None, title=None, LCI=None, UCI=None,org_prob=None):
    fig=plt.figure(figsize=(3,3))
    plt.title(title)

    #corner cases:
    data = [0 if x < 1e-2 else x for x in data]
    data = [1 if x > 0.99 else x for x in data]
    height, bins, patches = plt.hist(data, color='#49759C', edgecolor='#49759C', alpha=0.7)
    plt.fill_betweenx([0, height.max()], 0.01 if LCI<0.01 else LCI, 0.99 if UCI>0.99 else UCI, color='#49759C', alpha=0.2)  # Mark between 0 and the highest bar in the histogram
    plt.axvline(x=org_prob,color = 'k', linestyle = '-')
    plt.xlabel('HRD Probability')
    plt.ylabel('#replicates')
    plt.xticks([i/10 for i in range(11)])
    plt.xticks(rotation=90)
    fig=plt.tight_layout(pad=0.7)
    pp.savefig(plt.gcf())
