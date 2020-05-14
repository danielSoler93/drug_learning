import pandas as pd

def save_dict_to_df(dict):
    for name, vectors in dict.items():
        df = pd.DataFrame({"graph": vectors})
        df.to_csv(name, index=False)


