# Supplement to XXX et al. (in review)

## Contents

### Python Notebooks

 - `[00_NeptuneProcessing.ipynb](00_NeptuneProcessing.ipynb)`: Code to derive δ7Li values from raw Neptune data.
 - `01_DataCompilation.ipynb`: Code to compile isotopic, trace element and precipitation data, and unmix the composition of the experimental overgrowth from the seed crystals.
 - `02_Results.ipynb`: Plots and statistical analysis of the compiled data that are presented in the manuscript.

### Data

 - `data/`: contains all data files used in the analysis
    - `raw_Neptune_stripped.csv`: Output from Neptune analysis of Li isotope samples, reformatted from .xlsx to .csv with some header information removed.
    - `raw_conditions.csv`: Raw data containing precipitation conditions, trace element measurements and experimental metadata.
    - `raw_literature.csv`: A compilation of literature data containing paired Li/Ca and δ7Li analyses.
    - `processed_*.csv`, `processed_*.pkl`: Combined data files produced by `01_DataCompilation.ipynb` containing all solution and solid data. The `.pkl` files are in pickle format, and can be directly read by the `pandas` package, preserving all data structure and propagated uncertainty values.

### Helper Functions

 - `li_funks`: A collection of helper functions used in the notebooks to process and plot data.
