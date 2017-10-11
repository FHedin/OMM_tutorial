# Installation instructions for miniconda3

Create a conda environment using the provided YAML file:

```
conda env create -f environment.yml
```

Force install of the most recent version of nglview

```
conda upgrade -c bioconda --force nglview notebook ipywidgets 
```

Probably need to activate widgetsnbextension

```
python -m ipykernel install --sys-prefix
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter nbextension enable --py --sys-prefix nglview
```


