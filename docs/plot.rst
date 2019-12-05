Plot routines for a MCMax3D Model
*********************************

This module provides routines to plot the output of a MCMax3D model. 
This is just a starting point, things are likely to change dramaticaly. 


Usage example
-------------

The following example reads an MCMax3D model from the current working directory and makes a plot. 
(A matplotlib figure)

.. code-block:: python

   import mcmax3dpy.read as mread
   
   model=mread.read(modelDir=".",outDir="output")  

   fig=mplot.plot_cuts_zones(model.zones,"temp",rlim=25,zlim=[3,3500])


Source documentation
--------------------

.. automodule:: mcmax3dpy.plot

