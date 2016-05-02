#!/usr/bin/env python

# Author: Christian Vallentin <mail@vallentinsource.com>
# Website: http://vallentinsource.com
# Repository: https://github.com/MrVallentin/LinearMath
#
# Date Created: March 12, 2016
# Last Modified: May 02, 2016
#
# Developed and tested using Python 3.5.1

import unittest

from vec3 import *

class test_vec3(unittest.TestCase):
	def test_vec3(self):
		x, y, z = 0, 1, 2
		v = vec3(x, y, z)
		
		y = v.y = 4
		
		self.assertEquals(x, v.x)
		self.assertEquals(y, v.y)
		self.assertEquals(z, v.z)
	
	
	def test_vec3_to_list(self):
		# Turn vector into a list
		v = vec3(0, 1, 2)
		l = [ v[i] for i in range(len(v)) ]
		
		self.assertEquals(l, [ v.x, v.y, v.z ])
	
	
	def test_swizzling(self):
		v = vec3(0, 1, 2)
		
		self.assertEquals(v.xxx, vec3(v.x, v.x, v.x))
		self.assertEquals(v.zyx, vec3(v.z, v.y, v.x))
		
		self.assertEquals(v.rgb, vec3(v.x, v.y, v.z))
		self.assertEquals(v.rgb, vec3(v.r, v.g, v.b))
		self.assertEquals(v.xyz, vec3(v.r, v.g, v.b))
		
		self.assertEquals(v.xyzxyzxyz, (v.x, v.y, v.z, v.x, v.y, v.z, v.x, v.y, v.z))
	
	
	def test_compare(self):
		self.assertEquals(vec3(), vec3(0))
		self.assertEquals(vec3(), vec3(0, 0, 0))
		
		self.assertEquals(vec3(1), vec3(1, 1))
		self.assertEquals(vec3(1), vec3(1, 1, 1))
		
		self.assertNotEquals(vec3(1), vec3(1, 1, 2))
	
	
	def test_relational_operators(self):
		pass
	
	
	def test_subscript_operator(self):
		pass
	
	
	def test_arithmetic_operators(self):
		pass
	
