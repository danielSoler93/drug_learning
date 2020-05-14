from tqdm import tqdm
import pandas as pd
import drug_learning.datasets.CYP.provide_dataset as cyp
import drug_learning.two_dimension.Graph.graph as gp
import drug_learning.two_dimension.Graph.helpers as hp


def get_graph_from_CYP_dataset():
    #move to generators
    dataset = cyp.CYPDataset()
    train, test = dataset.get_X_Y()
    train_smiles = train["SMILES"].values
    test_smiles = test["SMILES"].values
    get_graph(train_smiles, test_smiles)

def get_graph(train_smiles, test_smiles):
    graph_train = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(train_smiles)]
    graph_test = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(test_smiles)]
    hp.save_dict_to_df({"graph_train.csv": graph_train, "graph_test.csv": graph_test})





if __name__ == "__main__":
    get_graph_from_CYP_dataset()
