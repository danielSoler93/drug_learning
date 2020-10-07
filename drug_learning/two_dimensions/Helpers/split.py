import os
#import glob
# from drug_learning.two_dimensions.Helpers import paralel as pl

def split_in_chunks(sdf_file, n_chunks):
    output_folder = f"split_sdfs_{n_chunks}"
    output_files = []
    print("Processing file", sdf_file)
    name, ext = os.path.splitext(sdf_file)
    name = name.replace("/", "_")
    mols = 0
    splits = 1
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    fw = open(os.path.join(output_folder, name+f"_{splits}"+ext), "w")
    output_files.append(os.path.join(output_folder, name+f"_{splits}"+ext))
    with open(sdf_file) as f:
        for line in f:
            fw.write(line)
            if line.startswith("$$$$"):
                mols += 1
            if mols == n_chunks:
                mols = 0
                fw.close()
                splits += 1
                fw = open(os.path.join(output_folder, name+f"_{splits}"+ext), "w")
                output_files.append(os.path.join(output_folder, name+f"_{splits}"+ext))
    return output_files

# if __name__ == "__main__":
#
#     files = glob.glob("ligands.sdf")
#     n_workers = 1
#     n_chunks = 1000
#     pl.paralelise_function(split_in_chunks, files, n_workers, n_chunks)
