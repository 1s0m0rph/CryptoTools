from crypto_tools.util_upper import *
from crypto_tools.ModularArith import *

"""
Rational maths
"""
from typing import Union,List

import numpy as np


class Rational:
	"""
	a rational number a/b where a and b are coprime* integers (b = 0 is thought of as being infinite)
	* iff rat_reduce is set to true
	"""

	FLOAT_APPROX_CONV_NUM = -1#what convergent will we use when converting from a float?

	'''
	if dec_display_digits is set to -1, we'll show as "a/b". otherwise we'll print as a decimal number out to that (max) precision
	'''
	def __init__(self,a:int,b:int,rat_reduce=True,dec_display_digits:int=10):
		if rat_reduce:
			g = ext_eucl_int(a,b,gcd_only=True)
			a //= g
			b //= g

		self.a = a
		self.b = b
		self.dec_display_digits = dec_display_digits

	def pairwise_check(self,other):
		if type(other) in [int,np.int32,np.int64]:
			other = Rational(other,1,rat_reduce=False,dec_display_digits=self.dec_display_digits)
		elif type(other) in [float,np.float64]:
			other = continued_frac_approx_convergents(other)[self.FLOAT_APPROX_CONV_NUM]
		return other

	def __add__(self,other):
		other = self.pairwise_check(other)
		return Rational(self.a*other.b+other.a*self.b,self.b*other.b,rat_reduce=True,
						dec_display_digits=self.dec_display_digits)

	def __neg__(self):
		return Rational(-self.a,self.b,dec_display_digits=self.dec_display_digits)

	def __sub__(self,other):
		return self+(-other)

	def __mul__(self,other):
		other = self.pairwise_check(other)
		return Rational(self.a*other.a,self.b*other.b,rat_reduce=True,dec_display_digits=self.dec_display_digits)

	def __abs__(self):
		return Rational(abs(self.a),abs(self.b),dec_display_digits=self.dec_display_digits)

	def __repr__(self):
		if self.dec_display_digits != -1:
			st = self.p_exp_st(self.dec_display_digits)
			return st
		else:
			if self.b == 0:
				return ("-" if self.a < 0 else "")+"INF"
			return "{}/{}".format(self.a,self.b)

	"""
	called with the unary "~" operator
	"""
	def __invert__(self):
		return Rational(self.b,self.a,dec_display_digits=self.dec_display_digits)

	def __truediv__(self,other):
		return self*(~other)

	def __floordiv__(self,other):
		return self/other

	def __le__(self,other):
		other = self.pairwise_check(other)
		diff = self-other
		return diff.a <= 0

	def __eq__(self,other):
		other = self.pairwise_check(other)
		return (self.a == other.a) and (self.b == other.b)

	def __lt__(self,other):
		return (self <= other) and (self != other)

	def __gt__(self,other):
		return not (self <= other)

	def __ge__(self,other):
		return not (self < other)

	def __pow__(self, power:int):#disallow noninteger powers so that we stay rational
		#rat reduction is unnecessary assuming we were already reduced
		#gcd(a,b) = 1 => gcd(a^p,b^p) = 1
		#because the prime fac. of a = p1^e1 * p2^e2 * ...
		#therefore a^p = (p1^e1 * p2^e2 * ...)^p = p1^(p*e1) * p2^(p*e2) * ...
		#and same for b^p, so if there were no common nonzero exponents beforehand there are still no common nonzero exponents
		if power < 0:
			return Rational(self.b**(-power),self.a**(-power),rat_reduce=False,
							dec_display_digits=self.dec_display_digits)  #(a/b)^-x = 1/(a/b)^x = 1/(a^x/b^x) = b^x/a^x
		if power == 0:
			return Rational(1,1,dec_display_digits=self.dec_display_digits)
		return Rational(self.a**power,self.b**power,rat_reduce=False,dec_display_digits=self.dec_display_digits)

	'''
	base-b expansion out to k "digits"
	output is a tuple: (integer part,list with noninteger coefficients for base 'base')
	'''
	def p_expansion(self,k,base=10):
		r = self
		i_part = self.a//self.b
		rats = []
		for i in range(k):
			c_prev = rats[i-1] if i > 0 else i_part
			r -= Rational(c_prev,base**i,dec_display_digits=self.dec_display_digits)
			bir = r * (base**(i+1))
			c_this = bir.a//bir.b
			rats.append(c_this)

		return i_part,rats

	def p_exp_st(self,k,base=10):
		ipart,nipart = self.p_expansion(k,base=base)
		#drop trailing zeroes
		while (len(nipart) > 0) and (nipart[-1] == 0):
			nipart = nipart[:-1]
		st = str(ipart)
		if len(nipart) > 0:
			if base <= 10:
				st += '.' + ''.join([str(x) for x in nipart])
			else:
				st += '.' + ''.join(['(' + str(x) + ')' for x in nipart])

		return st
	
	
	def cont_frac(self,final_elt_one=False):
		"""
		provide the continued fraction sequence for self:
		xvec = <x0,x1,x2,...,xk> such that
		x = x0 + 1/(x1 + 1/(x2 + 1/(... + 1/xk)))
		
		with the following properties:
			x0 is an integer
			xi (for all 0 < i <= k) is a positive integer
		
		if the optional final element =1 flag is set, xk is guaranteed to equal one (i.e. continue evaluating all the way). otherwise it is guaranteed to not equal one (unless k = 0)
		"""
		
		xvec = []
		# copy numerator/denominator since we're going to compute using these
		a = self.a
		b = self.b
		
		while a != 1:
			# LDA
			q,r = long_divide(a,b,non_neg_rem=True)
			# we now have q,r such that b*q + r = a
			# thus, the rat a/b can be written as q + (r/b) or q + 1/(b/r)
			# thus, q is the next element of the x-vector, and the next pair we'll perform the algorithm step on will be b/r
			xvec.append(q)
			a = b
			b = r
			#TODO I think we'd need one more step for when final elt one is true
		
		return xvec
	
	def disp_cf(self,final_elt_one=False):
		"""
		display self as a continued fraction
		"""
		
		raise NotImplementedError()
		#TODO
		
		pass

def st_idx_fuzzy(st,x):
	try:
		ret = st.index(x)
		return ret
	except ValueError:
		return None

'''
parse rational
'''
def rp(x:Union[str,int,float,Rational,List[int]],dec_display_digits=-1):
	if type(x) == int:
		return Rational(x,1,dec_display_digits=dec_display_digits)
	if type(x) == float:
		return Rational(1,1,dec_display_digits=dec_display_digits)*x
	if type(x) == Rational:
		return x
	if type(x) == str:
		decpt = st_idx_fuzzy(x,'.')
		if decpt is not None:
			return Rational(1,1,dec_display_digits=dec_display_digits)*float(x)
		divl = st_idx_fuzzy(x,'/')
		if divl is None:
			return Rational(int(x),1,dec_display_digits=dec_display_digits)

		#now do a proper parsing
		num = x[:divl]
		denom = x[divl+1:]
		return Rational(int(num),int(denom),dec_display_digits=dec_display_digits)
	if type(x) == list:
		# assume this is a cont frac list
		ratv = Rational(x[-1],1,dec_display_digits=dec_display_digits)
		for vec_elt in reversed(x[:-1]):
			# set equal to vec_elt + 1/<current value>
			ratv = Rational(1,1,dec_display_digits=dec_display_digits) / ratv
			ratv += Rational(vec_elt,1,dec_display_digits=dec_display_digits)
		
		return ratv

def eval_convergents(cidxs:List[int]):
	res = Rational(cidxs[-1],1)
	for cidx in cidxs[-2::-1]:
		res = ~res
		res += cidx
	return res

def continued_frac_convergents(r_inp:Rational) -> List[Rational]:
	#TODO is there a faster way than just brute forcing the actual convergents?
	r = abs(r_inp)
	i = r.a//r.b
	cidxs = [i]
	convs = [Rational(i,1)]
	rem = r - i
	while rem.a > 1:
		i = rem.b//rem.a
		rem = Rational(rem.b%rem.a,rem.a)
		cidxs.append(i)
		conv = eval_convergents(cidxs)
		convs.append(conv)
	convs.append(r)
	return convs


def continued_frac_approx_convergents(x:Union[float,np.float64],w=100) -> List[Rational]:
	if not np.isfinite(x):
		return [Rational(int(np.sign(x)),0)]
	#first generate a totally brain-dead guess (i.e. <integer part of x> + <rational part of x>*2^w / 2^w
	i = int(x)
	ratxnum = int((x-i)*(2**w))
	if ratxnum == 0:
		return [Rational(i,1)]
	rat = Rational(ratxnum,1<<w)+i
	convs = continued_frac_convergents(rat)
	return convs

