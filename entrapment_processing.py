"""Module for processing entrapment data"""
import re
import os
from collections import namedtuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pepxmls as pp


def clean_sequence(pep):
    """remove mod info and special chars from the peptide string"""
    clean = re.sub(f'{re.escape("[")}.*?{re.escape("]")}','', pep)
    clean = re.sub("-", "", clean)
    clean = re.sub('[0-9]', "", clean)
    return clean



def create_reference_dict(reference_csv):
    """generate reference dictionary for the ground truth files,
    spectrum: peptide"""
    dfr = pd.read_csv(reference_csv, header=None, names=["file", "peptide"])
    file_peptide = [(get_file_number(row[1]), clean_sequence(row[2])) for row in dfr.itertuples()]

    # generate dictionary
    reference_dict = {}
    for idx, peptide in file_peptide:
        reference_dict.setdefault(idx, []).append(peptide)

    return reference_dict


def get_file_number(item):
    """get file number from the dataframe"""
    return int(item.split('_')[-1])


def create_spectrast_pvalue_table(file_no, ref_dict):
    """Generate peptide table as an input to SpectraST
    when it constructs a spectral library with -pT option"""

    peptides = ref_dict[file_no]

    pep_formatted = []
    for charge in [2, 3, 4]:
        pep_formatted.extend([f"{peptide}/{charge}" for peptide in peptides])

    path = os.path.join("/home/dominik/entrapment/", f"pval_table_{file_no}.tsv")
    dfs = pd.DataFrame(pep_formatted)
    dfs.to_csv(path, sep='\t', header=False, index=False)


def evaluate_pvalues(negatives, to_evaluate, output_name, **kwargs):
    """generate p-values for entrapment PSMs based on
    the score distribution from decoy PSMs"""

    to_evaluate.sort()
    pvalue_iterator = (len([y for y in negatives if y >= x]) / negatives.size for x in to_evaluate)
    p_values = np.fromiter(pvalue_iterator, dtype=np.float32)

    with open(f'{output_name}.npy', 'wb') as file:
        np.save(file, p_values)
    #p_values = [get_empirical_pval(to_evaluate, x) for x in negs]

    if kwargs.get("plot"):
        plot_pvalues(p_values, kwargs.get("plot_option"))

    # return p_values


def get_empirical_pvalue(arr, x):
    """Find fraction of array that's larger than x"""

    length = next((i for i, element in enumerate(arr) if element >= x), len(arr))
    return (arr.size - length) / arr.size


def plot_pvalues(pvals, option):
    """plot the p-value distribution,
    ideally it should be a uniform distribution"""
    plt.rcParams.update({'font.size': 15})
    color = 'green' if 'decoy' not in option else 'orange'

    fig, _ = plt.subplots(figsize=(6,4))
    sns.histplot(pvals, stat='density', bins=25, color=color)
    plt.xlabel("p-value")
    plt.hlines(xmin=0, xmax=1, y=1, color='k', linestyles='dashed')
    fig.savefig(f"./p_values_{option}.png", dpi=600)


def plot_mass_distributions(data, tag_1, tag_2):
    """compare precursor mass distributions of
    negatives and entrapments"""

    fig, axes = plt.subplots(figsize=(5,5))

    for tag in [tag_1, tag_2]:
        for dec_label in [0, 1]:
            scores, emp_cdf = get_masses(data, tag, dec_label)
            plt.plot(scores, emp_cdf, label=f"{tag}_{dec_label}")

    axes.set_xlabel("precursor mass [Da]")
    axes.legend()
    axes.set_ylim(0, 1)

    fig.tight_layout()
    fig.savefig("./graphs/mass_distribution.png", dpi=500)


def get_masses(data, tag, dec_label):
    """get masses of rows with index containing
    the specified tag"""

    mass = data[data.index.str.contains(tag) & data.decoy == dec_label]['mass'].values
    length = len(mass)
    ran_vals = np.arange(1, length + 1)
    emp_cdf = ran_vals / length
    scores = sorted(mass)

    return scores, emp_cdf


def get_aa_mass():
	@@ -349,57 +156,268 @@ def get_aa_mass():

def peptide_mass(seq, aa_dict):
    """Calculate mass of a peptide"""
    mass = sum(aa_dict.get(aa, 0) for aa in seq)
    return round(mass, 2)


def pep_mass_hist(ref_dict, file_nos, aa_dict):
    """Plot mass histograms for list of files to be processed (based on pool.csv)"""

    plt.rcParams.update({'font.size': 13})
    fig, axs = plt.subplots(figsize=(5, 5))

    for file_no in file_nos:
        masses = [peptide_mass(pep, aa_dict) for pep in ref_dict[file_no]]
        sns.kdeplot(masses, label=file_no)

    axs.set_xlabel("peptide mass [Da]")
    axs.legend()
    fig.savefig("../graphs/pep_hist_all.png", dpi=500, bbox_inches="tight")


def plot_score_density(negs, decoys, ents):
    """Plots density plot of score distribution for
    negatives, decoys, entrapments"""
    # sns.set(rc={'figure.figsize': (5, 5)})

    plt.figure(figsize=(5, 5))
    sns.kdeplot(negs, label='negs')
    sns.kdeplot(decoys, label='decs')
    sns.kdeplot(ents, label='ents')

    plt.xlabel("TEV (Tide)")
    plt.legend()

    plt.savefig("../graphs/hist_testing.png", dpi=500, bbox_inches='tight')

    #sns.kdeplot(x[x.index.str.contains("3_") & x['label'] == True]['score'].to_numpy())
    #sns.histplot(x['score'], stat='density')



def get_pvalues_input(x, pos_tag, neg_tag, etag):
    """Returns score arrays for negatives, entrapments, decoys, positive samples"""

    negatives = x[(x.index.str.contains(f"{neg_tag}_")) &
             (~x['label']) &
             (x['decoy'] == 0)]['score'].to_numpy()

    entrapments = x[(x.index.str.contains(etag)) &
             (~x['decoy'])]['score'].to_numpy()

    decoys = x[(x['decoy'] == 1) &
               (~x['label']) &
               (x.index.str.contains(f"{pos_tag}_"))]['score'].to_numpy()

    positives = x[(x['label']) &
             (x.index.str.contains(f"{pos_tag}_"))]['score'].to_numpy()

    return negatives, entrapments, decoys, positives





class CompareResults:
    """class with methods to compare search results with decoy spectra
    vs. entrapment spectra"""

    def __init__(self, **kwargs):

        config_template = namedtuple("Configuration", kwargs.keys())
        self.config = config_template(*kwargs.values())
        self.ref_syn = create_reference_dict(self.config.ref_syn)
        self.results_df = None

    def run_analysis(self):
        """run the main analysis to get stats data
        from the entrapment vs decoy search results"""

        data = self.merge_pepxml_results()
        self.calculate_fdr_stats(data, mode=self.config.mode)

    def merge_pepxml_results(self):
        """read pepxml results, add ground truth labels
        and return combined dataframe"""

        positives_df, negatives_df, entrapments_df = self.read_pepxmls()

        entrapment_ratio = {"entrapment_ratio": len(entrapments_df)/len(positives_df[positives_df.protein == 1])}
        self.config.ratios.update(entrapment_ratio)

        positives_df = self.add_gt_labels(positives_df)
        negatives_df = self.add_gt_labels(negatives_df, 'neg')
        negatives_df = negatives_df[negatives_df['label'] != "remove"] # remove peptide shared with positives
        entrapments_df["label"] = len(entrapments_df) * (False,)

        all_d = self.combine_dataframes(positives_df, negatives_df, entrapments_df)
        all_d = self.add_main_score(all_d)
        all_d = self.add_counts(all_d)

        return all_d


    def add_main_score(self, df):
        """Add main score to the dataframe"""

        df['score'] = df[self.config.score]

        if self.config.engine == 'SpectraST' and self.config.score == 'p-value':
            df['score'] = -0.02 * np.log(df['p_value'] * df['hits_num'] / 500)

        return df


    def calculate_fdr_stats(self, data, mode='PSM'):
        """calculating FDP for FDR thresholds
        provided by decoys and entrapments"""
        data.sort_values("score", ascending=False, inplace=True)
        data = self.add_counts(data)

        if mode == "peptide":
            idx = data.groupby(['peptide', 'entrapment'])['score'].transform(max) == data['score']
            data = data[idx]

        data['targets'] = data['targets'].cumsum()

        for fdr_type in ('decoy', 'entrapment', 'GT'):
            data = self.add_fdr_columns(data, fdr_type)

        self.results_df = data
        return data


    def add_fdr_columns(self, data, opt='decoy'):
        """add column with FDR estimates to the main dataframe"""

        if opt == 'GT':
            data[f'{opt}_fdr'] = data[opt].cumsum()/data['targets']
            return data

        if opt == 'decoy':
            data[f'{opt}_fdr'] = data[opt].cumsum()/data['targets']
            decoy_filter = (data.decoy == 1) & (data.entrapment == 0)
            non_decoy_filter = (data.decoy == 0) & (data.entrapment == 0)
            decoy_ratio = len(data[decoy_filter]) / len(data[non_decoy_filter])
            self.config.ratios.update({f"{opt}_ratio": decoy_ratio})

        if opt == 'entrapment':
            entrap_filter = (data[opt] == 1) & (data.decoy == 0)
            data[f'{opt}_fdr'] = data[entrap_filter][opt].cumsum() / data['targets']

        data[f'{opt}_fdr'] = data[f'{opt}_fdr']/self.config.ratios[f"{opt}_ratio"]
        return data


    def add_counts(self,dfr):
        """add counts of decoys, entraps, and incorrect targets
        to be used in FDR/FDP calculation"""

        dfr["decoy"] = dfr.apply(lambda row: 1 if row['protein'] == 0
                                 and self.config.etag not in row[0]
                                 else 0, axis=1)
        dfr['entrapment'] = dfr.apply(lambda row: 1
                                      if self.config.etag in row[0]
                                      else 0, axis=1)
        dfr['GT'] = dfr.apply(lambda row: 0
                              if row['label']
                              else 1, axis=1)
        dfr['targets'] = dfr.apply(lambda row: 1
                                   if row['protein'] == 1
                                   else 0, axis=1)

        return dfr

    def read_pepxmls(self):
        """Reads content of all input pepXML files"""
        return [pp.PepXMLResults(x).read_pepxml()
                for x in self.config.input_files]


    def _get_gt_labels(self, df_row, option=''):
        """check whether peptide found in synthetic data is the correct one"""
        bare_seq = clean_sequence(df_row[1])

        if option == 'neg':
            shared_peps = set(self.ref_syn[self.config.pos_tag]).intersection(set(self.ref_syn[self.config.neg_tag]))
            return "remove" if bare_seq in shared_peps else False

        if self.config.pos_tag in self.ref_syn:
            return bare_seq in self.ref_syn[self.config.pos_tag]

        return False


    def add_gt_labels(self, dfr, option=""):
        """add new columns with ground truth labels to the dataframe"""
        if option == 'neg':
            dfr['label'] = [self._get_gt_labels(row, 'neg') for row in dfr.itertuples()]
        else:
            dfr['label'] = [self._get_gt_labels(row) for row in dfr.itertuples()]
        return dfr


    @staticmethod
    def combine_dataframes(*dataframes):
        """combine dfs"""
        return pd.concat(dataframes, axis=0, ignore_index=False)


    def add_fdp_decoy_fdr(self, dfs):
        """Calculates FDP and decoy-based FDR estimates"""
        new_df = dfs[~dfs.index.str.contains(self.config.etag)].copy()
        df_arange = np.arange(1, len(new_df) + 1)
        new_df.loc[:,'GT_fdr'] = new_df['GT'].cumsum() / df_arange
        new_df.loc[:,'decoy_fdr'] = 2 * new_df['decoy'].cumsum() / df_arange
        return new_df


    def add_entrapment_fdr(self, original_df):
        """Calculates entrapment-based FDR estimate"""
        new_df = original_df[original_df.decoy == 0].copy()
        new_df['ents_fdr'] = self.config.ratios["entrapment_ratio"] * new_df['entrapment'].cumsum() / np.arange(1, len(new_df) + 1)
        return new_df

    @staticmethod
    def plot_fdp_vs_fdr_results(decoy_df, entrapment_df, output_name, font_size):
        """Plots FDP vs FDR estimation results"""

        plt.rcParams.update({'font.size': font_size})
        fig, axs = plt.subplots(figsize=(5,5))

        axs.plot(decoy_df['decoy_fdr'], decoy_df['GT_fdr'], label="decoy")
        axs.plot(entrapment_df['ents_fdr'], entrapment_df['GT_fdr'], label="entrapment")
        axs.plot([0, 0.05], [0, 0.05], color='grey')

        axs.set_xlim(0, 0.05)
        axs.set_ylim(0, 0.05)
        axs.set_xlabel("FDR")
        axs.set_ylabel("FDP")
        axs.legend()

        fig.savefig(f"../graphs/{output_name}.png", dpi=600, bbox_inches="tight")


    def plot_fdp_fdr_results(self, dfs, suffix):
        """Execute processing and plotting of FDP vs FDR results"""
        #ent_ratio = 1/(len(x[x.index.str.contains(en_tag)])/len(x[~x.index.str.contains(en_tag)]))
        dec_df = self.add_fdp_decoy_fdr(dfs)
        ent_df = self.add_entrapment_fdr(dfs)

        outname = f"pos_26_neg_42_ent_{self.config.etag}_{self.config.mode}_{suffix}"

        self.plot_fdp_vs_fdr_results(dec_df, ent_df, outname, font_size=15)


    def evaluate_plot_pvalues(self, outname):
        """evaluate and plot p-value distributions for entrapments and decoys"""

        negatives, entrapments, decoys, _ = get_pvalues_input(self.results_df,
                                                            self.config.pos_tag,
                                                            self.config.neg_tag,
                                                            self.config.etag)
        evaluate_pvalues(negatives, decoys, f'{outname}_decoy_pvalues',
                                 plot=True, plot_option=f'{outname}_decoy')
        evaluate_pvalues(negatives, entrapments, f'{outname}_entrapment_pvalues',
                                 plot=True, plot_option=f'{outname}_entrapment')