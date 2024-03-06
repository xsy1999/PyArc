# PyArc
A python package of method for computing radiative recombination coefficients and absorption coefficients from first principles based on scheme proposed by [Xie Zhang *et al.*](https://doi.org/10.1021/acsenergylett.8b01297). The code also offers a serious of functions to output the calculation results graphically and analyze the mechanism beneath those results. More details about those functions can be found in our recent paper. 

## Installation
NONRAD is implemented in python and can be installed through `pip`.
Dependencies are kept to a minimum and include standard packages such as `numpy`, `scipy`, and `matplotlib`.

#### With pip
As always with python, it is highly recommended to use a virtual environment.
Install directly from github,
```
$ pip install git+https://github.com/xsy1999/PyArc
```

#### For development
To install NONRAD for development purposes, clone the repository
```
$ git clone https://github.com/xsy1999/PyArc
```
then install the package in editable mode with development dependencies
```
$ pip install -r requirements.txt
```

## Usage
A tutorial notebook will be added later
The basic steps are summarized below:

0. Perform a first-principles calculation of the target superconductor system. The transmatrix values (Transmatrix file) can be obtain with setting the parameter 'LOPTICS = True' based on a inplemented VASP version. A good explanation of its methodology can be found in this [Recombination in Semiconductors](https://doi.org/10.1017/CBO9780511470769). A high quality calculation with enough k sampling points in Brillouin Zone when start electronic self-consistency steps is necessary as input for the radiative recombination coefficients since too coarse original k-grid sampling in Brillouin zone will introduce some deviations.
1. Calculate the potential energy surfaces for the configuration coordinate diagram. This is facilitated using the `get_cc_structures` function. Extract the relevant parameters from the configuration coordinate diagram, aided by `get_dQ`, `get_PES_from_vaspruns`, and `get_omega_from_PES`.
2. Calculate the electron-phonon coupling matrix elements, using the method of your choice (see [our paper]() for details on this calculation with `VASP`). Extraction of the matrix elements are facilitated by the `get_Wif_from_wavecars` or the `get_Wif_from_WSWQ` function.
3. Calculate scaling coefficients using `sommerfeld_parameter` and/or `charged_supercell_scaling`.
4. Perform the calculation of the nonradiative capture coefficient using `get_C`.


#### How to Cite
If you use our code to calculate nonradiative capture rates, please consider citing
```
@article{Xie Zhang_first-principles_2018,
	title = {First-Principles Analysis of Radiative Recombination in Lead-Halide Perovskites},
	volume = {3},
	doi = {10.1021/acsenergylett.8b01297},
	number = {10},
	journal = {ACS Energy Lett.},
	author = {Alkauskas, Audrius and Yan, Qimin and Van de Walle, Chris G.},
	month = sep,
	year = {2018},
	pages = {2329},
}
```
