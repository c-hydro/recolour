Recolour
===========================
Welcome to the **Recolour** GitHub repository. The **REprocess paCkage for sOiL mOistUre pRoducts** 
is used for preventing and reducing the hydrogeological risk.


Background
**********


Components
**********


Prerequisites
*************
In order to use the Flood-PROOFS forecasting chain, users are strongly raccomanted to control if the following characteristics, libraries and packages are available and correctly installed on their machine.

Usually, Flood-PROOFS libraries are installed on **Linux Debian/Ubuntu 64bit** environment and all libraries, packages and applications must be compilled and/or installed in according with this operative system.

All codes, subroutines and scripts are developed using both **Python** (version 3 and greater) [2_] and **Fortran** (version 2003 and greater) [3_]. QGIS geographic information system (version 2.18 and greater) [4_] is used to develop tools to organize, create and control static and dynamic datasets. R (version 3.4.4 and greater) [5_] is used to perform statistical analysis.

The libraries and the packages are mainly divided in four categories:

    • python3 packages and applications;
    • R packages and applications;
    • fortran libraries;
    • other software and applications (Jupyter Notebook, QGIS, Panoply, cdo, ncview ...).

The section for installing all needed libraries and environments is usually named **fp-envs** and the users can find it in Flood-PROOFS
modelling system repository hosted by GitHub [1_].

Python3 libraries
-----------------

The python3 standard library is not sufficient to correctly install all Flood-PROOFS applications; for this reason some extra libraries are needed to guarantee all functionalities. 
To install all python3 libraries a bash script named **setup_fp_env_python.sh** is provided [6_]; basically, the script calls a **miniconda** [7_] installation that allow to get all needed libraries and install them into “$HOME/user/fp_libs_python/” folder. During the installation, a virtual environment named “fp_env_python” is created too.
Once all libraries are correctly installed and configurated, to activate “fp_env_python” by command-line is necessary to execute the following:

.. code-block:: bash
    
   >> export PATH="$HOME/fp_libs_python/bin:$PATH"
   >> source activate fp_env_python

By default, the **fp_env_python** environment is shown in parentheses () or brackets [] at the beginning of your command prompt:

.. code-block:: bash

   (fp_env_python) >> 

Activating the virtual enviroment permits to use a correct configuration andall applications and scripts of Flood-PROOFS forecasting chain will work properly.


Potential Users
***************
The Flood-PROOFS Modelling System has been released to enable different applications (for example local/regional scenario assessment) and further development by external users.

Potential users are anticipated to predominately be interested in the ability to run the system with local data (including scenario modelling) and to modify the system with new capabilities. The potential collaborators have expressed a range of potential goals for their use of the modelling system, including performing comparisons with existing models, tailoring the hydrological performance to specific land uses and cropping types.

Broadly speaking, there are four potential user categories of the FloodPROOFS modelling system:

    • **Data user**: who accessing the model outputs for using them in their analysis.
    • **Case study user**: who work to evaluate his/her case using data over a selected time period.
    • **Applying users**: who would primarily be interested in applying the current model to a region of interest using localised and/or scenario data where available.
    • **Contributor users**: who will extend the capabilities of the model with new research and coding (modify the system with new capabilities)

It is expected that the majority of early adopters of the FloodPROOFS modelling system will be Applying users looking to apply the system with local data/scenarios, with more Contributor users adopting the system as it becomes well known and established.


Contribute and Guidelines
*************************
We are happy if you want to contribute. Please raise an issue explaining what is missing or if you find a bug. We will also gladly accept pull requests against our master branch for new features or bug fixes.

If you want to contribute please follow these steps:

    • fork the one of the Flood-PROOFS repositories to your account;
    • clone the repository, make sure you use "git clone --recursive" to also get the test data repository;
    • make a new feature branch from the repository master branch;
    • add your feature;
    • please include tests for your contributions in one of the test directories;
    • submit a pull request to our master branch.

Authors
*******
All authors involved in the library development for the Flood-PROOFS modelling system are reported in this authors_ file.


License
*******
By accessing or using the Flood-PROOFS modelling system, code, data or documentation, you agree to be bound by the FloodPROOFS license available. See the license_ for details. 


Changelog
*********
All notable changes and bugs fixing to this project will be documented in this changelog_ file.


References
**********
| [1_] CIMA Hydrology and Hydraulics GitHub Repository
| [2_] Python programming language
| [3_] QGIS project
| [4_] Conda environment manager
| [5_] Hydrological Model Continuum codes

.. _1: https://github.com/c-hydro
.. _2: https://www.python.org/
.. _3: https://github.com/c-hydro/fp-env
.. _4: https://conda.io/miniconda.html
.. _5: https://github.com/c-hydro/hmc-dev
.. _license: LICENSE.rst
.. _changelog: CHANGELOG.rst
.. _authors: AUTHORS.rst
