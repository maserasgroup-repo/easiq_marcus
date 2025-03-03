# Easiq_marcus
Python package that facilitates the design of the Gaussian calculations and the curation of data to obtain energy barriers for SET/EnT process. 

## Getting Started

### Prerequisities 
* python >=3.6

*  `pip` or `conda` (for package management)

### Download and Installing the Code
Easiq_marcus is a collection of multiple Python scripts. The recommended steps to download and install the package are:

1. **Clone the repository**:
  ```bash
   git clone https://github.com/maserasgroup-repo/easiq_marcus.git
   cd easiq_marcus
  ```
2.  **Set up the virtual environment** \
2.1.  It is recommended to create the environment with the `.yml` file because all the dependencies are installed. \
Using conda (it can also be set up with pip and venv):
  ```bash
   conda env create -f easiqm.yml
   conda activate easiqm
   pip install -e .
  ```

2.2  Alternative, manually onstall dependencies:
If you prefer not to use the `.yml` file, you can manually create a Conda environment and install the dependencies:

  ```bash
   conda env create -name easiqm python=3.9
   conda activate easiqm
   pip install -e .
  ```

3. **Install additional required libraries:** \
If the step 2.2 is carried out, you will need to manually install the following dependencies:
\
The easiq_marcus package relies on the `pyssian` and `pyssianutils` libraries for managing Gaussian input and output files. These libraries must be installed before using easiq_marcus:

  ```bash
  pip install git+https://github.com/maserasgroup-repo/pyssian.git@develop
  pip install git+https://github.com/maserasgroup-repo/pyssian-utils.git
  ```

The remaining prerequisites are Python's built-in modules and do not require separate installation.

4. **Install the easiq_marcus package:**\
If you have already followed 1 + 2.1 or 2.2, the package is already installed, and you can skip this step. Otherwise, install the package directly from the GitHub repository:

  ```bash
pip install -e git+https://github.com/maserasgroup-repo/easiq_marcus.git
  ```

***(!)For private repository access:*** \
If the repository is private, organization members can install the package using their GitHub username and token:

  ```bash
pip install -e git+https://<USERNAME>:<TOKEN>@github.com/maserasgroup-repo/easiq_marcus.git#egg=easiq_marcus
  ```

5.  **If you want to use the Jupyter notebook, install the IPython kernel:** \
Dependencies for running the Jupyter notebook:
- jupyter
- ipykernel

```bash
pip install jupyter ipykernel
python -m ipykernel install --user --name=easiqm --display-name "ieasiqm"
```
Type jupyter-notebook to access to the files. \
Unzip the folder `calculations.zip` and move the files to the folder `tutorials_easiq_marcus` to use `4point_method.ipynb` and `harsphere_method.ipynb`


## easiq_marcus package

There are four python scripts (contained in the PyOSET directory). The functionality of each script is detailed above.  Additionally, it is worth noting that the code is thoroughly documented and its functionalities are explained in the aforementioned article **(url article)**.

* **utilities.py**: script that provides the tools to submit calculations and read the dictionaries. 

* **get_info.py**: script that extracts the information from the Gaussian output files of the reactants and products and save it in a dictionary. This dictionary is used in the next script.

* **4points.py**: script that creates input files required to obtain the value of the free energy of the SET/EnT activation barrier using the approximation of the 4-points from two parabolas.

* **get_results.py**: script that provides the free energy of the SET/EnT transfer following approximations of the Marcus theory.

## Example Usage - Obtain the energetic parameters involved in the SET/EnT step



### Authors
Lucía Morán-González - [lmoranglez](https://github.com/lmoranglez) \
Albert Solé-Daura - [asoledaura](https://github.com/asoledaura) \
Feliu Maseras - [maserasgroup](https://github.com/maserasgroup)

### Citation
**Add article citation**


