from rdkit import Chem
from rdkit.Chem import Draw
from sklearn import decomposition
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from PIL import ImageFont
from PIL import ImageDraw

df = pd.read_csv("schrodinger_fragment_collection_epik_MorganFP.csv", index_col=0)

pca = decomposition.PCA(n_components=2)
pca.fit(df)
pca_res = pca.transform(df)
x = pca_res[:,0]
y = pca_res[:,1]

mols = Chem.SDMolSupplier("schrodinger_fragment_collection_epik.sd")

frag_imgs = list(map(Chem.Draw.MolToImage, mols))
mol_names = [mol.GetProp("_Name") for mol in mols]

for i, img in enumerate(frag_imgs):
    draw = ImageDraw.Draw(img)
    font = ImageFont.truetype("Roboto-Black.ttf", 17)
    draw.text((50, 250), mol_names[i] ,"black",font=font)

img_array = list(map(np.asarray, frag_imgs))

# create figure and plot scatter
fig = plt.figure()
ax = fig.add_subplot()
line, = ax.plot(x,y, ls="", marker="o")

imagebox = OffsetImage(img_array[0], zoom=0.5)
imagebox.image.axes = ax

xybox=(50., 50.)
ab = AnnotationBbox(imagebox, (0,0), xybox=xybox, xycoords='data',
        boxcoords="offset points",  pad=0.1,  arrowprops=dict(arrowstyle="-"))
# add it to the axes and make it invisible
ax.add_artist(ab)
ab.set_visible(False)

def hover(event):
    # if the mouse is over the scatter points
    if line.contains(event)[0]:
        # find out the index within the array from the event
        ind = line.contains(event)[1]["ind"]
        # get the figure size
        w,h = fig.get_size_inches()*fig.dpi
        ws = (event.x > w/2.)*-1 + (event.x <= w/2.)
        hs = (event.y > h/2.)*-1 + (event.y <= h/2.)
        #if event occurs in the top or right quadrant of the figure,
        #change the annotation box position relative to mouse.
        ab.xybox = (xybox[0]*ws, xybox[1]*hs)
        # make annotation box visible
        ab.set_visible(True)
        # place it at the position of the hovered scatter point
        ab.xy =(x[ind[0]], y[ind[0]])
        # set the image corresponding to that point
        imagebox.set_data(img_array[ind[0]])
    else:
        #if the mouse is not over a scatter point
        ab.set_visible(False)
    fig.canvas.draw_idle()

# add callback for mouse moves
fig.canvas.mpl_connect('motion_notify_event', hover)
plt.show()
