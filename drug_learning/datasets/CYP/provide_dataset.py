import pandas as pd
import os

DIR = os.path.dirname(os.path.abspath(__file__))
TRAIN = os.path.join(DIR, "CYP2C9_dataset_training.csv")
TEST = os.path.join(DIR, "CYP2C9_dataset_testing.csv")
DATASET = [TRAIN, TEST]

class CYPDataset:

    def __init__(self):
        self.train, self.test = DATASET

    def get_X_Y(self):
        """
        Return training and test set for CYPs
        :return:
            train-pandas_dataframe
            test-pandas_dataframe
        """
        return  pd.read_csv(self.train), pd.read_csv(self.test)