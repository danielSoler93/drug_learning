

class Saver():
    def __init__(self, df):
        self.df = df

    def to_csv(self, output):
        self.df.to_csv(output)

    def to_parquet(self, output):
        self.df.to_parquet(output, compression='gzip')

    def to_feather(self, output):
        self.df = self.df.reset_index() #Requiered by pandas. Parquet is a better option
        self.df.to_feather(output)

    def to_hdf(self, output):
        self.df.to_hdf(output, key='df', mode='w')

    def to_pickle(self, output):
        self.df.to_pickle(output)