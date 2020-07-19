# drug_learning
Open package with DL and ML tools for drug discovery. It is based in a 3 branch package: 2D-QSAR, 3D-QSAR, simulation+AI

## To developers

master: This should always be working and ready for production!

devel: Only for pre-releases. This should be the test branch before merging to master. **All custom function should have their own unit-test here!**

2D: ligand-based, 2D Deep learning methodologies to screen billions of compounds (Joan)

3D: structural-based,  3D Deep learning methodologies to retrieve patterns out of the output of simulations (Alexis)

predictions: python-API to predict over commercial libraries of compounds (Ana)

simulation: Scripts to analyse the output of the simulation (sub-pocket, hbonds, pharmacophore...) (Carles&Marti)

To avoid merge conflicts on the maste branch, when merging please push to your branch  and then pull request to devel  and to  master **(your_branch > devel > master)** 


