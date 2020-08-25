From python API
=================

Transform molecules into fingerprints
--------------------------------------
Create an object of a specific fingerprint class. Then use the ``fit`` method to load the sdf file, and ``transform``
to obtain the fingerprints.

::

  from drug_learning.two_dimensions.Input import fingerprints as fp

  input_sdf = "../datasets/ligands.sdf"

  mor = fp.MorganFP()
  structures = mor.fit(input_sdf)
  features = mor.transform()


Set the applicability domain
-------------------------------------

::

  from drug_learning.two_dimensions.Input import applicability_domain as ad

  AD = ad.ApplicabilityDomain()
  AD.fit(train_features)
  AD.predict(test_features)
