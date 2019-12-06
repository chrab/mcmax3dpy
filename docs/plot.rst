Plot routines for a MCMax3D Model
*********************************

This module provides routines to plot the output of a MCMax3D model. 
This is just a starting point, things are likely to change significantly in the future. 


Usage example
-------------

The following example reads an MCMax3D model from the current working directory and makes a plot. 
(more precisely a matplotlib figure). The easiest way to test this is a Jupyter notebook

.. code-block:: python

   import mcmax3dpy.read as mread
   
   model=mread.read(modelDir=".",outDir="output")  

   fig=mplot.plot_cuts_zones(model.zones,"temp",rlim=25,vlim=[3,3500])


Source documentation
--------------------

.. automodule:: mcmax3dpy.plot

