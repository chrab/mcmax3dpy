Read routines for a MCMax3D Model
*********************************

This module provides several routines to read the output of a MCMax3D model. 
All the data belonging to a MCMax3D model is put into an hierachical data structure (:class:`~mcmax3dpy.read.DataMCMax3D`).

The module provides a routine to read "all" the data of a MCMax3D model and spezialized routines  
to read only distinct model data (e.g. only the SED of a model).

Usage example
-------------

The following example reads the zones files of a MCMax3D model from the current working directory.
There are also more spezialed read routines (see :mod:`mcmax3dpy.read`). 

.. code-block:: python

   import mcmax3dpy.read as mread
   import mcmax3dpy.plot as mplot
   
   model=mread.read(modelDir=".",outDir="output")  
   print(model)
   

Source documentation
--------------------

.. automodule:: mcmax3dpy.read

