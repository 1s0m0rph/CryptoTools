"""
This file defines the common interface into the project, as well as a few generified functions that apply in more than one place
"""

from util_upper import *
from classic_crypto import *
from FiniteFields import *
from ModularArith import *
from factoring import *
from RSA import *
from EllipticCurves import *

def ext_eucl(a,b,gcd_only=False):
	if type(a) == FiniteFieldPoly:
		return FFP_ext_eucl(a, b, gcd_only=gcd_only)
	elif type(a) in [int,ModInteger]:
		return ext_eucl_int(a,b,gcd_only=gcd_only)
	else:
		raise AttributeError("Unkown extended-euclid type: {}".format(type(a)))

def gcd(a,b):
	return ext_eucl(a,b,gcd_only=True)