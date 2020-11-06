import numpy as np
from scipy.spatial import distance


class ApplicabilityDomain():

    def __init__(self, tresh=5):
        self.x_train = None
        self.x_test = None
        self.fitted = None
        self.tresh = tresh

    def fit(self, x_train):
        self.x_train = x_train
        distances = np.array([distance.cdist([x], self.x_train) for x in self.x_train])
        distances_sorted = [np.sort(d[0]) for d in distances]
        d_no_ii = [ d[1:] for d in distances_sorted]
        k = int(round(pow(len(self.x_train), 1/3)))

        d_means = [np.mean(d[:k][0]) for d in d_no_ii] #medium values
        Q1 = np.quantile(d_means, .25)
        Q3 = np.quantile(d_means, .75)
        IQR = Q3 - Q1
        d_ref = Q3 + 1.5*(Q3-Q1) #setting the reference value
        n_allowed =  []
        all_allowed = []
        for i in d_no_ii:
            d_allowed = [d for d in i if d <= d_ref]
            all_allowed.append(d_allowed)
            n_allowed.append(len(d_allowed))

        #selecting minimum value not 0:
        min_val = [np.sort(n_allowed)[i] for i in range(len(n_allowed)) if np.sort(n_allowed)[i] != 0]

        #replacing 0's with the min val
        n_allowed = [n if n!= 0 else min_val[0] for n in n_allowed]
        all_d = [sum(all_allowed[i]) for i, d in enumerate(d_no_ii)]
        self.thresholds = np.divide(all_d, n_allowed) #threshold computation
        self.thresholds[np.isinf(self.thresholds)] = min(self.thresholds) #setting to the minimum value where infinity
        self.fitted = True
        return self.thresholds

    def predict(self, x_test):
        print(f"Molecules before applicability domain filtering {x_test.shape[0]}")
        assert isinstance(x_test, np.ndarray), "x_test needs to be a numpy array"
        assert self.fitted, "Need to perform object.fit(x_train) before"
        self.x_test = x_test
        test_names= ["sample_{}".format(i) for i in range(self.x_test.shape[0])]
        d_train_test = np.array([distance.cdist([x], self.x_train) for x in self.x_test])
        self.n_insiders = []

        for i, name in zip(d_train_test, test_names): # for each sample
            idxs = [j for j,d in enumerate(i[0]) if d <= self.thresholds[j]] #saving indexes of training with threshold < distance
            self.n_insiders.append(len(idxs))

        self.n_insiders = np.array(self.n_insiders)
        idx = np.array([True if insiders > self.tresh else False for insiders in self.n_insiders])
        inside_domain = x_test[idx, :]
        print(f"Molecules after applicability domain filtering {inside_domain.shape[0]}")
        return inside_domain

    def save(self, output="ad"):
        np.save(output, self)

    def load(self, object):
        return np.load(object, allow_pickle=True).item()


