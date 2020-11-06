import joblib

class Model():

    def __init__(self, file):
        self.file = file
        self.obj = self.load()

    def load(self):
        loaded_model = joblib.load(self.file)
        return loaded_model