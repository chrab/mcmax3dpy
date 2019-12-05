# README #

This package is still under development. 
At some point it should be become a usefull tool to read in and plot output of the 3D dust radiative transfer code MCMax3D (developed by Michiel Min). 
Any contributions are welcome!

Requirements
============
mcmax3dpy uses several additional python packages which are commonly used in the astronomical community. 
If you use [anaconda](https://www.anaconda.com/distribution/) all this packages should be available in your python distribution. 
The following packages are required

* *matplotib* required for the plotting part only, version>=2 is recommended  
* *astropy*     version >=2 is recommended
* *numpy*       no known special requirements
* *scipy*       no known special requirements

If you use the setup script (see Installation) those packages will be installed automatically if 
they are not included in your python distribution. We stronglhy recommend to use python3.

Installation
============
Currently the easiest way to use it is to clone this repository and install the package directly from the source:

* change into a directory of your choice and 
* clone the repository 


```
  git clone https://github.com/chrab/mcmax3dpy.git
``` 
 
* change into the newly created mcmax3dpy directory and type:


```
python setup.py develop
```

This will install the package in your current python environment. 
The develop options allows to update the python code (e.g. via git) without the need to reinstall the package.

If you do not have root access to install python packages, this should work

```
python setup.py develop --user
```

Code Update
===========
Simply type 

```
git pull 
```

in the mcmax3dpy directory. You can directly use the updated code (no reinstall required).

Documentation
=============
Please check out the documentation! Click on the badge!

[![Documentation Status](https://readthedocs.org/projects/mcmax3dpy/badge/?version=latest)](https://mcmax3dpy.readthedocs.io/en/latest/?badge=latest)
