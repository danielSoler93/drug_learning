from tqdm import tqdm
import pandas as pd
import drug_learning.datasets.CYP.provide_dataset as cyp
import drug_learning.two_dimension.Graph.graph as gp


def run():
    dataset = cyp.CYPDataset()
    train, test = dataset.get()
    train_smiles = train["SMILES"].values
    test_smiles = test["SMILES"].values
    molecule_vectors_train = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(train_smiles)]
    molecule_vectors_test = [gp.GraphBuilder(smile).build_graph() for smile in tqdm(test_smiles)]
    for vectors, name in zip([molecule_vectors_train, molecule_vectors_test], ["graph_train.csv", "grsaph_test.csv"]):
        df = pd.DataFrame({"graph": vectors})
        df.to_csv(name, index=False)






if __name__ == "__main__":
    run()
