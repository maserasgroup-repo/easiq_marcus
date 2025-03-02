# Easiq_marcus
Python package that facilitates the design of the Gaussian calculations and the curation of data to obtain energy barriers for SET/EnT process. 

## Getting Started

### Prerequisities 
* python >=3.6

*  `pip` or `conda` (for package management)

* easiq_marcus package manages Gaussian input and output files using the library **pyssian** and **pyssianutils** (available at https://github.com/maserasgroup-repo/pyssian and https://github.com/maserasgroup-repo/pyssian-utils, respectively). The libraries must be installed prior to using PyOSET package for manipulating Gaussian files. It is recommended to follow the specified installation procedure to ensure the correct branch installation. Previous versions require to install the develop branch. \

* The remaining prerequisities are Python's built-in modules.

### Download and Installing the Code
Easiq_marcus is a collection of multiple Python scripts. The recommended steps to follow are:

1. **Clone the repository**:
  ```bash
   git clone https://github.com/maserasgroup-repo/easiq_marcus.git
   cd easiq_marcus
  ```
2.  **It is recommended to set up a virtual environment with the .yml file**.
  Using conda (it can also be set up with pip and venv):
  ```bash
   conda env create -f easiqm.yml
   conda activate easiqm
  ```

3.  **Alternatively to the point (2), it can be installed as bellow:**
  ```bash
   conda env create -name easiqm
   conda activate easiqm
  ```
4   **Install the repository**
  ```bash
  pip install -e git+https://github.com/maserasgroup/easiq_marcus.git
  ```
  while the repository is private, organization people can install the package using 
  $ 
  ```bash
   pip install -e git+https://<USERNAME>:<TOKEN>@github.com/maserasgroup-repo/easiq_marcus.git#egg=easiq_marcus
  ```
4.  **If you want to use the Jupyter notebook, install the IPython kernel:** \
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


