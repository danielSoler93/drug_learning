import drug_learning.datasets.CYP.provide_dataset as cyp
import drug_learning.two_dimension.Graph.graph as gp

def run():
    train, test = cyp.CYPDataset().get()
    smiles = train["SMILES"].values
    graph = gp.GraphBuilder(smiles[0])
    graph.build_graph()




if __name__ == "__main__":
    run()
