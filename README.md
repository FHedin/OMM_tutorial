# Installation instructions for miniconda3

Create a conda environment using the provided YAML file:

```
conda env create -f environment.yml
```
If you get an error message for nglview, it is probably needed to activate some extensions:

```
python -m ipykernel install --sys-prefix
jupyter nbextension enable --py --sys-prefix widgetsnbextension
jupyter nbextension enable --py --sys-prefix nglview
```

