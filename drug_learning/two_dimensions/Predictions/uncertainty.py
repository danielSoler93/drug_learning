import numpy as np

class InterModel():

    def __init__(self, model, fingerprints):
        self.model = model
        self.fingerprints = fingerprints

    def get_uncertainty_intermodel(self):
        pred1 = self.model.named_estimators_.svc1.predict(self.fingerprints)
        pred2 = self.model.named_estimators_.svc2.predict(self.fingerprints)
        pred3 = self.model.named_estimators_.svc3.predict(self.fingerprints)

        scores = []
        for p1, p2, p3 in  zip(pred1, pred2, pred3):
            results = [p1, p2, p3]
            zeros = len([r for r in results if r == 0])
            ones = len([r for r in results if r == 1])
            score = max([zeros, ones])
            scores.append(score)
        return scores

class IntraModel():

    def __init__(self, preds):
        self.preds = preds

    def get_uncertainty_intramodel(self):
        return [[0, 1]]*len(self.preds)


class Uncertainty(InterModel, IntraModel):

     def __init__(self, preds, model, fingerprints):
         InterModel.__init__(self, model, fingerprints)
         IntraModel.__init__(self, preds)


     def filter_by_uncertanty(self, tresh_intra=0.7, tresh_inter=3):
         print(f"Molecules before uncertainty filter {self.preds.shape[0]}")
         intra = self.get_uncertainty_intramodel()
         inter = self.get_uncertainty_intermodel()
         filtered = []
         for i, (pred, score, p) in enumerate(zip(self.preds, inter, intra)):
             if score >= tresh_inter and max(p)>tresh_intra:
                 filtered.append(0)
             else:
                 filtered.append(1)
         filtered = np.array(filtered)
         molecules_after_filter = self.preds[np.where(filtered == 0)]
         print(f"Molecules after uncertainty filter {molecules_after_filter.shape[0]}")
         return  molecules_after_filter
