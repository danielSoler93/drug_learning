import setuptools

setuptools.setup(name='drug_learning',
      version="1.0.0",
      url = "https://github.com/danielSoler93/drug_learning",
      description='Convert 3D molecules into fingerprints or descriptors.',
      maintainer='Ana Sanchez',
      maintainer_email='ana.sanchez01@estudiant.upf.edu',
      packages=setuptools.find_packages(),
      install_requires=["numpy", "pandas", "scipy", "tables", "pyarrow", "fastparquet", "mordred", "pytest"],
      classifiers=[
       "Programming Language :: Python :: 3",
       "License :: OSI Approved :: MIT License",
       "Operating System :: OS Independent",
       ],
     )
