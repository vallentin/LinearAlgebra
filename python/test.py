#!/usr/bin/env python

# Author: Christian Vallentin <mail@vallentinsource.com>
# Website: http://vallentinsource.com
# Repository: https://github.com/MrVallentin/LinearMath
#
# Date Created: March 12, 2016
# Last Modified: May 02, 2016
#
# Developed and tested using Python 3.5.1

import os, sys

import unittest
import glob

#from linearmath import *

if __name__ == "__main__":
	test_cases = glob.glob("./linearmath/test_*.py".replace("/", os.sep))
	test_cases = [ test_case[:-3].replace(os.sep, ".").lstrip(".") for test_case in test_cases ]
	
	print("Test Cases:\n- " + "\n- ".join(test_cases))
	print(" ")
	
	sys.path.append("./linearmath".replace("/", os.sep))
	
	suite = unittest.defaultTestLoader.loadTestsFromNames(test_cases)
	
	#result = unittest.TextTestRunner().run(suite)
	result = unittest.TextTestRunner(verbosity=2).run(suite)
