"""
This file defines the common interface into the project, as well as a few generified functions that apply in more than one place
"""

from crypto_tools.util_upper import *
from crypto_tools.classic_crypto import *
from crypto_tools.FiniteFields import *
from crypto_tools.ModularArith import *
from crypto_tools.factoring import *
from crypto_tools.RSA import *
from crypto_tools.Isogenies import *
from crypto_tools.EllipticCurves import *
from crypto_tools.Lattices import *
from crypto_tools.RationalArith import *

def ext_eucl(a,b,gcd_only=False):
	if type(a) == FiniteFieldPoly:
		return FFP_ext_eucl(a, b, gcd_only=gcd_only)
	elif type(a) in [int,ModInteger]:
		return ext_eucl_int(a,b,gcd_only=gcd_only)
	else:
		raise AttributeError("Unkown extended-euclid type: {}".format(type(a)))

def gcd(a,b):
	return ext_eucl(a,b,gcd_only=True)