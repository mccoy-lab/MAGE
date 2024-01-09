## GeoVar Figure Generation

The scripts in this directory generate the following figures in the main manuscript:

* 


### Running Figure Generation

To create the appropriate figures, we recommend the following steps:


1. Create the conda environment:

```
conda env create -f env.yaml
conda activate mage-geovar
```

2. Run the primary scripts for table generation: 

```
python3 finemapping_metaCS_dapg_kgpex_replicate.py
```

NOTE: this step specifically creates our meta credible sets to allow comparison between GTEx multi-tissue results and MAGE results.  

3. Create figures 

```
python3 finemapping_metaCS_geovar_plots.py
```

