# PyArc
A python package of method for computing radiative recombination coefficients and absorption coefficients from first principles based on scheme proposed by [Xie Zhang *et al.*](https://doi.org/10.1021/acsenergylett.8b01297). The code also offers a series of functions to output the calculation results graphically and analyze the mechanism beneath those results. More details about those functions can be found in our recent paper. 

## Installation
PyArc is implemented in python and can be installed through `pip`.
Dependencies are kept to a minimum and just include standard packages such as `numpy`, `scipy`, and `matplotlib`.

#### With pip
As always with python, it is highly recommended to use a virtual environment.
Install directly from github,
```
$ pip install git+https://github.com/xsy1999/PyArc
```

#### For development
To install PyArc for development purposes, clone the repository
```
$ git clone https://github.com/xsy1999/PyArc
```
then install the package in editable mode with development dependencies
```
$ pip install -r requirements.txt
```

## Usage
A tutorial notebook will be added later.
The basic steps are summarized below:

0. Perform a first-principles calculation of the target superconductor system. The transmatrix values (Transmatrix file) can be obtain with setting the parameter 'LOPTICS = True' using a modified VASP version. A good explanation of its methodology can be found in this book [Recombination in Semiconductors](https://doi.org/10.1017/CBO9780511470769). A high quality calculation with enough k sampling points in Brillouin Zone when start electronic self-consistency steps is necessary for the radiative recombination coefficients since too coarse original k-grid sampling in Brillouin zone will introduce some deviations.
1. Calculate a denser k-grid of both eigenvalues and transmatrix values by the interpolation method offered in our code. This is facilitated using the `Interp_main.py` function. Through setting relevant parameters like the magnification factor or valence bands and conduction bands, the interpolated eigenvalues and transmatrix values will be obtained as `Eigen_geninterp.dat` and `matrix_fine.dat` file respectively.
 The way to excute `Interp_main.py` function can be
 ```
 $python Interp_main.py -VB 34 35 36 -CB 37 38 39 -m 10 10 10
 ```
2. Calculate the absorption coefficients and radiative recombination coefficients. This is facilitated by the `Coefficients_main.py` function. Also, for both two coefficients calculations, a series of parameters should be delivered to the function according to demands as below. Then files like `Absorption-GaAs.dat` or `Radiative-GaAs_trail.dat` will be generated to store those results.
 ```
 $python Coefficients_main.py -A True -R True -E 1 3 300 # absorption coefficients
 $python Coefficients_main.py -R True -u 3.89 -V 180 -S false -T 100 200 300 # radiative recombinaiton coefficients
 or 
 $python Coefficients_main.py -A True R True -u 1.58 -V 180 -S false -T 100 200 300 # calculate both
 ```

3. Plot the calculation results including absorption coefficients and radiative recombination coefficients under different conditions. This is facilitated by the `Figure_main.py` function. It also support to plot carrier density distribution. The way to execute this function can be
 ```
 $python Figure_main.py -A True -R True -VB 34 35 36 -CB 37 38 39 -E 1.42
 ```
 users can obtain final results files under the work path they defined

## Contributing
To contribute, see the above section on installing [for development](#for-development).
Contributions are welcome and any potential change or improvement should be submitted as a pull request on [Github](https://github.com/xsy1999/PyArc/).
Potential contribution areas are:
 - [ ] implement a command line interface
 - [ ] add more robust tests for various functions

#### How to Cite
If you use our code to calculate absorption coefficients or radiative recombination coefficients, please consider citing
```
@article{Xie Zhang_first-principles_2018,
	title = {First-Principles Analysis of Radiative Recombination in Lead-Halide Perovskites},
	volume = {3},
	doi = {10.1021/acsenergylett.8b01297},
	number = {10},
	journal = {ACS Energy Lett.},
	author = {Xie Zhang, Jimmy-Xuan Shen, Wennie Wang, and Chris G. Van de Walle},
	month = sep,
	year = {2018},
	pages = {2329},
}
```
