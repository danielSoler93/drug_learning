import pandas as pd

class Fingerprint():

    def __init__(self, file):
        self.file = file
        self.format = self.get_format()
        self.df = self.load()

    def load(self):
        if self.format == "csv":
            return pd.read_csv(self.file)
        elif self.format == "pq":
            return pd.read_parquet(self.file)
        elif self.format == "pkl":
            return pd.read_pickle(self.file)

    def get_format(self):
        return self.file.rsplit(".")[-1]