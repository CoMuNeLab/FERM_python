# FERM_python
Simulate the Feature-Enriched Radiation Model (FERM) on Python.
Code includes parallelization of the main routine.

## Installation
Open a terminal and run
```bash
python -m pip install .
```
or, if you want to install the package locally (develop mode), run
```bash
python -m pip install --editable .
```

## Getting started
See `./tests/` for examples of usage

- *sigma*: standard deviation of the gaussian distribution to sample the absorbance/absorption threshold
- *nb_particules*: number of particules to produce for each pair of points
- *path_niche_array*: matrix of the value the human suitability for each point
- *path_x*: longitude coordinates
- *path_y*: latitude coordinates
- *path_pop*: .tif file where each geospatial point is associated with a population value
- *save_mobility*: where to save the mobiltiy matrix

## License
[MIT](https://choosealicense.com/licenses/mit/)
