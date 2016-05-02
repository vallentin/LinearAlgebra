#!/usr/bin/env python

# Author: Christian Vallentin <mail@vallentinsource.com>
# Website: http://vallentinsource.com
# Repository: https://github.com/MrVallentin/LinearMath
#
# Date Created: March 12, 2016
# Last Modified: May 02, 2016
#
# Developed and tested using Python 3.5.1

from os.path import dirname, basename, isfile
import glob

modules = glob.glob(dirname(__file__) + "/*.py")
modules = [ basename(f)[:-3] for f in modules if isfile(f) ]

__all__ = [ f for f in modules if not f.startswith("test_") and not f.endswith("__init__") ]

del modules
