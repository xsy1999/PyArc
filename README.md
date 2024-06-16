# PyArc
This is a python package for computing absorption and radiative coefficients from first principles based on the scheme proposed by [Xie Zhang *et al.*](https://doi.org/10.1021/acsenergylett.8b01297). The code also offers a series of functions to visualize the calculation results graphically and analyze the microscopic mechanism beneath those results. More details about those functions can be found in our recent paper. 

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
A Jupyter notebook is available in the main directory to demonstrate the use of the code for a specific example of GaAs.

The basic steps are summarized below:

0. Perform a first-principles calculation for the target material system. The transmatrix values (Transmatrix file) can be obtained by setting the parameter 'LOPTICS = True' in the INCAR file while using a modified VASP version. The concrete patch file for VASP can be found in the vasp_patch directory and we offered file for different versions.

    To apply the patch file, please download the related 'optics.diff' file into VASP installation directory and use the patch command
    ```
    $patch ./src/linear_optics.F < optics.diff
    ```
    and then recompile VASP. After those steps, one can find the "Transmatrix" file in the calculation directory.

    A good explanation of our methodology can be found in this book [Recombination in Semiconductors](https://doi.org/10.1017/CBO9780511470769). A high-quality first-principles calculation with enough k points for sampling the Brillouin Zone in the electronic self-consistent loop is necessary for accurate absorption and radiative coefficients since too coarse original k-grid in the Brillouin zone may introduce pronounced deviations.

1. Calculate the eigenvalues and transmatrix for a dense k-point grid by using the interpolation method offered in our code. This is facilitated using the `Interp_main.py` function. Through setting relevant parameters like the magnification factor and the valence and conduction band indices to interpolate, the interpolated eigenvalues and transmatrix values will be obtained as `Eigen_geninterp.dat` and `matrix_fine.dat` file respectively.
    The way to excute `Interp_main.py` function is as follows (default parameters for testfiles):
    ```
    $python Interp_main.py
    ```
    and also could add some parameters like
    ```
    $python Interp_main.py -VB 34 35 36 -CB 37 38 39 -m 10 10 10
    ```
2. Calculate the absorption and radiative  coefficients. This is facilitated by the `Coefficients_main.py` function. Also, for the calculations of both coefficients, a series of parameters should be delivered to the function according to demands as below. Then files like `Absorption-GaAs.dat` or `Radiative-GaAs_trail.dat` will be generated to store those results.
run with default parameters for test files
    ```
    $python Coefficients_main.py
    ```
    or add parameters like
    ```
    $python Coefficients_main.py -A True -R True -E 1 3 300 # absorption coefficients
    $python Coefficients_main.py -R True -u 3.89 -V 180 -S false -T 100 200 300 # radiative recombinaiton coefficients
    or 
    $python Coefficients_main.py -A True R True -u 1.58 -V 180 -S false -T 100 200 300 # calculate both
    ```

3. Plot the calculation results for absorption and radiative coefficients under different conditions. This is facilitated by the `Figure_main.py` function. It also supports visualization of the carrier density distribution. The way to execute this function is as follows (with default parameters for testfiles):
    ```
    $python Figure_main.py
    ```
    or set parameters like
    ```
    $python Figure_main.py -A True -R True -VB 34 35 36 -CB 37 38 39 -E 1.42
    ```
    users can obtain final results files under the work path they defined
 
4. Users can find the GaAs data as an example in the directory of 'test_files'. And all input files such as KPOINTS can be found in 'input_files' under the 'test_files' directory. The Transmatrix file and Wannier interpolated data files are zipped. The size of the zipped Wannier file is large and please download it seperately.

## Contributing
To contribute, see the above section on installing [for development](#for-development).
Contributions are welcome and any potential change or improvement should be submitted as a pull request on [Github](https://github.com/xsy1999/PyArc/).
Potential contribution areas are:
 - [ ] implement a command line interface
 - [ ] add more robust tests for various functions

#### How to Cite
If you use our code to calculate absorption or radiative  coefficients, please consider citing
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
