"""Processing of pepXML files"""
from pyteomics import pepxml
import pandas as pd



class PepXMLResults:
    """Use synthetic spectra to verify whether entrapment queries can be used as a proxy for
real false positives, then use decoy counting to verify"""
    def __init__(self, filename) -> None:
        self.spectra = {}
        self.filename = filename
        self.score_names = []


    def read_pepxml(self):
        """read pepxml file and return
        pandas dataframe with some PSM features"""
        data = pepxml.read(self.filename)

        filtered = [x for x in data if "search_hit" in x.keys()]
        top_only = [x for x in filtered if len(x["search_hit"]) == 1]

        self.spectra = {}
        [self.get_spec_items(item) for item in top_only]

        sel_data = pd.DataFrame.from_dict(self.spectra, orient="index",
                                      columns=["peptide", *self.score_names, "protein", "mass"])

        return sel_data

    def get_spec_items(self, spec):
        """retrieve some useful data from each spectrum
        search results"""
        top_hit = spec["search_hit"][0]

        self.score_names = top_hit["search_score"].keys()
        score_vals = top_hit["search_score"].values()

        mass = spec["precursor_neutral_mass"]
        peptide = top_hit["peptide"]
        protein = 1
        # if a decoy protein, assign 0 instead
        if str.lower(top_hit["proteins"][0]["protein"][0]) == "d":
            protein = 0
        value = f"{peptide}"
        key = spec["spectrum"]
        self.spectra[key] = (value, *score_vals, protein, mass)