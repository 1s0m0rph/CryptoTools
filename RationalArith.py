from util_upper import *
from ModularArith import *



class Rational:
	"""
	a rational number a/b where a and b are coprime* integers (b = 0 is thought of as being infinite)

	* iff rat_reduce is set to true
	"""

	def __init__(self,a: int,b: int,rat_reduce=True):
		if rat_reduce:
			g = ext_eucl_int(a,b,gcd_only=True)
			a //= g
			b //= g

		self.numerator = a
		self.b = b

	def pairwise_check(self,other):
		if type(other) in [int,np.int32,np.int64]:
			other = Rational(other,1,rat_reduce=False)
		elif type(other) in [float,np.float64]:
			other = continued_frac_approx_convergents(other)[-1]
		return other

	def __add__(self,other):
		other = self.pairwise_check(other)
		return Rational(self.numerator*other.b+other.numerator*self.b,self.b*other.b,rat_reduce=True)

	def __neg__(self):
		return Rational(-self.numerator,self.b)

	def __sub__(self,other):
		return self+(-other)

	def __mul__(self,other):
		other = self.pairwise_check(other)
		return Rational(self.numerator*other.numerator,self.b*other.b,rat_reduce=True)

	def __abs__(self):
		return Rational(abs(self.numerator),abs(self.b))

	def __repr__(self):
		if self.b == 0:
			return ("-" if self.numerator < 0 else "")+"INF"
		return "{}/{}".format(self.numerator,self.b)

	"""
	called with the unary "~" operator
	"""

	def __invert__(self):
		return Rational(self.b,self.numerator)

	def __truediv__(self,other):
		return self*(~other)

	def __floordiv__(self,other):
		return self/other

	def __le__(self,other):
		other = self.pairwise_check(other)
		diff = self-other
		return diff.numerator <= 0

	def __eq__(self,other):
		other = self.pairwise_check(other)
		return (self.numerator == other.numerator) and (self.b == other.b)

	def __lt__(self,other):
		return (self <= other) and (self != other)

	def __gt__(self,other):
		return not (self <= other)

	def __ge__(self,other):
		return not (self < other)

def eval_convergents(cidxs:List[int]):
	res = Rational(cidxs[-1],1)
	for cidx in cidxs[-2::-1]:
		res = ~res
		res += cidx
	return res

def continued_frac_convergents(r_inp:Rational) -> List[Rational]:
	#TODO is there a faster way than just brute forcing the actual convergents?
	r = abs(r_inp)
	i = r.numerator//r.b
	cidxs = [i]
	convs = [Rational(i,1)]
	rem = r - i
	while rem.numerator > 1:
		i = rem.b//rem.numerator
		rem = Rational(rem.b%rem.numerator,rem.numerator)
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
	rat = Rational(ratxnum,1<<w) + i
	convs = continued_frac_convergents(rat)
	return convs