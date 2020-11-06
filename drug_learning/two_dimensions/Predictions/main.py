import numpy as np
import argparse
from drug_learning.two_dimensions.Input.applicability_domain import ApplicabilityDomain
from drug_learning.two_dimensions.Predictions.uncertainty import Uncertainty
from drug_learning.two_dimensions.Predictions.model import Model
from drug_learning.two_dimensions.Predictions.fingerprints import Fingerprint



class Predictor(ApplicabilityDomain):

    def __init__(self, model_file, fingerprint_file, ad_file, ad_tresh=5, intra_tresh=0.7, inter_tresh=3):
        self.ad_tresh = ad_tresh
        self.intra_tresh = intra_tresh
        self.inter_tresh = inter_tresh
        self._model = Model(model_file)
        self._fingerprints = Fingerprint(fingerprint_file)
        self._applicability_domain = ApplicabilityDomain().load(ad_file)

    def predict(self):
        inside_domain = self._applicability_domain.predict(self._fingerprints.df.values)
        predictions = self._model.obj.predict(inside_domain)
        self.uncertaunty = Uncertainty(predictions, self._model.obj, self._fingerprints.df)
        molecules_after_filter = self.uncertaunty.filter_by_uncertanty(self.intra_tresh, self.inter_tresh)
        self._report(molecules_after_filter)
        return molecules_after_filter

    def _report(self, molecules):
        actives = np.where(molecules == 1)[0].shape[0]
        #inactives = np.where(molecules == 0)[0].shape[
        print(f"Molecules predicted as actives: {actives}")

def parse_args(parser):
    parser.add_argument('model', type=str,
                        help='model file')
    parser.add_argument('fingerprints', type=str,
                        help='fingerprints to predict file')
    parser.add_argument('applicability_domain', type=str,
                        help='applicability domain file')
    return parser



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser = parse_args(parser)
    args = parser.parse_args()
    pred = Predictor(args.model, args.fingerprints, args.applicability_domain,
                     inter_tresh=3)
    molecules = pred.predict()


