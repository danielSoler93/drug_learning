from tqdm import tqdm
import pandas as pd
import drug_learning.datasets.CYP.provide_dataset as cyp
import drug_learning.two_dimension.Graph.graph as gp
import drug_learning.two_dimension.Graph.helpers as hp


def run():
    #move to generators
    dataset = cyp.CYPDataset()
    train, test = dataset.get()
    train_smiles = train["SMILES"].values
    test_smiles = test["SMILES"].values
    molecule_vectors_train = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(train_smiles)]
    molecule_vectors_test = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(test_smiles)]
    hp.save_dict_to_df({"graph_train.csv": molecule_vectors_train, "graph_test.csv": molecule_vectors_test})





if __name__ == "__main__":
    run()
