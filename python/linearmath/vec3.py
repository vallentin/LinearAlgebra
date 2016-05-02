#!/usr/bin/env python

# Author: Christian Vallentin <mail@vallentinsource.com>
# Website: http://vallentinsource.com
# Repository: https://github.com/MrVallentin/LinearMath
#
# Date Created: March 12, 2016
# Last Modified: May 02, 2016
#
# Developed and tested using Python 3.5.1

import math
import operator

#from vec2 import vec2

class vec3(object):
	x, y, z = 0, 0, 0
	
	def __init__(self, x = 0, y = None, z = None):
		if isinstance(x, vec3):
			self.x += x.x
			self.y += x.y
			self.z += x.z
		#elif isinstance(x, vec2):
		#	self.x += x.x
		#	self.y += x.y
		elif hasattr(x, "__getitem__"):
			self.x += x[0]
			self.y += x[1]
			self.z += x[2]
		else:
			if y == None:
				y = x
			if z == None:
				z = y
			
			self.x += x
			self.y += y
			self.z += z
	
	
	def __str__(self):
		return "vec3(%.2f, %.2f, %.2f)" % (self.x, self.y, self.z)
	
	def __repr__(self):
		return "vec3(%r, %r, %r)" % (self.x, self.y, self.z)
	
	
	
	
	def _toindex(self, key):
		key = str(key).lower()
		if key == "0" or key == "x" or key == "r" or key == "s":
			return 0
		elif key == "1" or key == "y" or key == "g" or key == "t":
			return 1
		elif key == "2" or key == "z" or key == "b" or key == "p":
			return 2
		else:
			return -1
	
	def _toname(self, index):
		if index == 0:
			return "x"
		elif index == 1:
			return "y"
		elif index == 2:
			return "z"
		else:
			return None
	
	
	def __len__(self):
		return 3
	
	def __setitem__(self, key, value):
		index = self._toindex(key)
		if index == 0:
			self.x = value
		elif index == 1:
			self.y = value
		elif index == 2:
			self.z = value
		else:
			raise IndexError("Invalid subscript '" + str(key) + "' to " + type(self).__name__)
	
	def __getitem__(self, key):
		index = self._toindex(key)
		if index == 0:
			return self.x
		elif index == 1:
			return self.y
		elif index == 2:
			return self.z
		else:
			raise IndexError("Invalid subscript '" + str(key) + "' to " + type(self).__name__)
	
	
	def __setattr__(self, name, value):
		for letter in name:
			attr_name = self._toname(self._toindex(letter))
			
			if attr_name is None:
				raise AttributeError(type(self).__name__ + " has no attribute '" + letter + "'")
			
			super().__setattr__(attr_name, value)
	
	def __getattr__(self, name):
		name = name.lower()
		name_len = len(name)
		
		swizzling = []
		for letter in name:
			swizzling.append(self[letter])
		
		if name_len == 1:
			return swizzling[0]
		#elif name_len == 2:
		#	return vec2(swizzling)
		elif name_len == 3:
			return vec3(swizzling)
		#elif name_len == 4:
		#	return vec4(swizzling)
		else:
			#return swizzling
			return tuple(swizzling)
	
	
	def __eq__(self, other):
		if hasattr(other, "__getitem__") and len(other) == 3:
			return self.x == other[0] and self.y == other[1] and self.z == other[2]
		else:
			return False
	
	def __ne__(self, other):
		if hasattr(other, "__getitem__") and len(other) == 3:
			return self.x != other[0] or self.y != other[1] or self.z != other[2]
		else:
			return True
	
	
	def _op(self, other, f):
		if isinstance(other, vec3):
			return vec3(f(self.x, other.x), f(self.y, other.y), f(self.z, other.z))
		elif hasattr(other, "__getitem__"):
			return vec3(f(self.x, other[0]), f(self.y, other[1], f(self.z, other[2])))
		else:
			return vec3(f(self.x, other), f(self.y, other), f(self.z, other))
	
	def _rop(self, other, f):
		if hasattr(other, "__getitem__"):
			return vec2(f(other[0], self.x), f(other[1], self.y), f(other[2], self.z))
		else:
			return vec2(f(other, self.x), f(other, self.y), f(other, self.z))
	
	def _iop(self, other, f):
		if hasattr(other, "__getitem__"):
			self.x = f(self.x, other[0])
			self.y = f(self.y, other[1])
			self.z = f(self.z, other[2])
		else:
			self.x = f(self.x, other)
			self.y = f(self.y, other)
			self.z = f(self.z, other)
		return self
	
	
	def __add__(self, other):
		return self._op(other, operator.add)
	def __radd__(self, other):
		return self._rop(other, operator.add)
	def __iadd__(self, other):
		return self._iop(other, operator.add)
	
	
	def __sub__(self, other):
		return self._op(other, operator.sub)
	def __rsub__(self, other):
		return self._rop(other, operator.sub)
	def __isub__(self, other):
		return self._iop(other, operator.sub)
	
	
	def __mul__(self, other):
		return self._op(other, operator.mul)
	def __rmul__(self, other):
		return self._rop(other, operator.mul)
	def __imul__(self, other):
		return self._iop(other, operator.mul)
	
	
	def __div__(self, other):
		return self._op(other, operator.div)
	def __rdiv__(self, other):
		return self._rop(other, operator.div)
	def __idiv__(self, other):
		return self._iop(other, operator.div)
	
	
	def __floordiv__(self, other):
		return self._op(other, operator.floordiv)
	def __rfloordiv__(self, other):
		return self._rop(other, operator.floordiv)
	def __ifloordiv__(self, other):
		return self._iop(other, operator.floordiv)
	
	
	def __truediv__(self, other):
		return self._op(other, operator.truediv)
	def __rtruediv__(self, other):
		return self._rop(other, operator.truediv)
	def __itruediv__(self, other):
		return self._iop(other, operator.truediv)
	
	
	def __mod__(self, other):
		return self._op(other, operator.mod)
	def __rmod__(self, other):
		return self._rop(other, operator.mod)
	
	
	def __divmod__(self, other):
		return self._op(other, operator.divmod)
	def __rdivmod__(self, other):
		return self._rop(other, operator.divmod)
	
	
	def __pow__(self, other):
		return self._op(other, operator.pow)
	def __rpow__(self, other):
		return self._rop(other, operator.pow)
	
	
	def __neg__(self):
		return vec3(operator.neg(self.x), operator.neg(self.y), operator.neg(self.z))
	
	def __pos__(self):
		return vec3(operator.pos(self.x), operator.pos(self.y), operator.pos(self.z))
	
	def __abs__(self):
		return vec3(operator.abs(self.x), operator.abs(self.y), operator.abs(self.z))
	
	def __invert__(self):
		return vec3(-self.x, -self.y, -self.z)
	
	
	def clone(self):
		return vec3(self)
	
