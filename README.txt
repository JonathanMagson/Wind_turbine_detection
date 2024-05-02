Wind turbine a contrario detection on Sentinel-2 images
==============================================================

Version 0.1 - September 16, 2021
by nicolas mandroux <nimandroux@gmail.com>


Introduction
------------

This is an implementation of an algorithm for wind turbine detection on Sentinel-2 images.
optical satellite time series.  This code is part of the following publication
and was subject to peer review:

  "[[COMPLETER]"


Files
-----

README.txt    			- This file
requirements.txt    		- Required python modules 
main_IPOL_OH_nogdal.py 	- Main python function, used for execution
utils_eol_compute.py  		- Auxiliary python function
utils_eol_meta_nogdal.py	- Auxiliary python function
run_miguel.sh          	- Executable for IPOL demo
test.tif			- Test image
test_meta.txt			- Test image's metadata
input_0.png			- Conversion of test.tif in .png
nfa.tif			- Nfa map produced by algorithm applied on test.tif
nfa.png			- Nfa map converted in .png
nfa_detec.png			- Detection map produced by algorithm


Running
-------

A typical execution is:

  python3 main_IPOL_OH_nogdal.py --filename test.tif --metaname test_meta.txt --t_NFA 1 --t_shadow 25 --t_hub 50

The input images can only be in TIFF format.


Copyright and License
---------------------

Copyright (c) 2021 nicolas mandroux <nimandroux@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.



