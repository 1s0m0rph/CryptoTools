from FiniteFields import *

class EllipticPoint:

	"""
	A point on the elliptic curve specified by the given polynomial
	"""

	'''
	Takes in a finite field polynomial as parameter, but thinks of it as being y^2 = <that>
	for_factor is a parameter which tells us that this EC point is intended to be used to factor a number
		this changes doubling behavior to not throw an exception on impossible division
		(generally this means p will not be prime)
	'''
	def __init__(self,poly:FiniteFieldPoly,x,y,for_factor=False,factor_canary=False):
		assert(poly.dgr == 3)
		self.is_inf = False
		self.for_factor = for_factor
		self.factor_canary = factor_canary#this is used for when we're factoring so we can detect (and return) non-coprime with p items

		if (not factor_canary) and ((x is None) or (y is None)):#we'll take this to mean infinity
			self.is_inf = True
		else:
			self.x = ModInteger(x,poly.p)#for factoring, this will contain the item
			self.y = ModInteger(y,poly.p)#as will this
		self.poly = poly
		self.p = poly.p
		if (not factor_canary) and (not self.is_inf) and (self.y**2 != poly[self.x]):
			raise AttributeError("Elliptic point {} is not on the curve specified by y² = {}".format((x,y),poly))


	def __repr__(self):
		if self.is_inf:
			return 'INF'
		else:
			return 'E{}: ({}, {})'.format(self.p,self.x.x,self.y.x)#simple -- could make more complex later


	def pairwise_check(self,other):
		if type(other) != EllipticPoint:
			return EllipticPoint(self.poly,*other,for_factor=self.for_factor)#think of it is a list/tuple
		if other.poly != self.poly:
			raise AttributeError("Comparing elliptic points on different curves")
		return other

	def __eq__(self, other):
		other = self.pairwise_check(other)
		if self.is_inf:
			return other.is_inf
		return (self.x == other.x) and (self.y == other.y)

	def __add__(self, other):
		other = self.pairwise_check(other)
		if self.is_inf:
			return other
		if other.is_inf:
			return self
		if self == other:
			#tangent line special case
			#E: y^2 = ax^3 + bx^2 + cx + d
			#dE: 2yy' = 3ax^2 + 2bx + c
			#dE: dy/dx = (3ax^2 + 2bx + c)/(2y) [plug in for x and y to get slope]
			if self.y == 0:
				return EllipticPoint(self.poly,None,None,for_factor=self.for_factor)#infinity

			denominator = self.y*2
			eval_poly = FiniteFieldPoly(self.p,[self.poly.coef[0]*3,self.poly.coef[1]*2,self.poly.coef[2]])
			if self.for_factor:
				try:
					m = eval_poly[self.x]/denominator
				except ZeroDivisionError:
					#this indicates we have found something which has nontrivial gcd with the modulus p
					#give us that thing
					return EllipticPoint(self.poly,denominator,denominator,for_factor=True,factor_canary=True)
			else:
				m = eval_poly[self.x]/denominator#if this throws an error anyway the user fukked up
			#point-slope for equation
			#y - y1 = m(x - x1)
			#y = m(x - x1) + y1
			#equate this ^^ to our poly
			#(mx - mx1 + y1)^2 = ax^3 + bx^2 + cx + d
			#m^2x^2 - 2m^2x1 + 2mxy1 + m^2x1^2 - 2mx1y1 + y1^2 = ax^3 + bx^2 + cx + d
			#only really care about the x^2 term anyway which is (m^2 - b)x^2
			#this means that xr = -(m^2 + b - 2x1)
			xr = -(self.x*2 - m**2 + self.poly.coef[1])
			#and yr is given by the above formoola
			yr = -(m*(xr - self.x) + self.y)#don't forget to negate

			return EllipticPoint(self.poly,xr,yr,for_factor=self.for_factor)
		elif self.x == other.x:
			#"a line passes through inf iff it is vertical"
			return EllipticPoint(self.poly,None,None,for_factor=self.for_factor)#infinity

		#the line going through both of us is given by
		denominator = other.x - self.x
		if self.for_factor:#bad division can also happen here
			try:
				m = (other.y-self.y)/denominator
			except ZeroDivisionError:
				return EllipticPoint(self.poly, denominator, denominator, for_factor=True, factor_canary=True)
		else:
			m = (other.y-self.y)/denominator
		#(y - y1) = m(x - x1)
		#y = mx - mx1 + y1
		#c = (y1 - mx1)
		c = self.y - m*self.x

		#this intersects the curve at (mx+c)^2 = ax^3 + bx^2 + cx + d (as well as us and the other point)
		#m^2x^2 + 2mcx + c^2 = ax^3 + bx^2 + cx + d
		#-(b - m^2) = sum of roots = x1 + x2 + xr
		#xr = -(x1 + x2 + b - m^2)
		xr = -(self.x + other.x + self.poly.coef[1] - m**2)
		#use the y = mx + c eqn to find y
		yr = m*xr + c
		#actual y is negated
		yres = -yr
		return EllipticPoint(self.poly,xr,yres,for_factor=self.for_factor)


	def __neg__(self):
		if self.is_inf:
			return self #-inf = inf (inf is additive identity)

		#otherwise negate just y
		return EllipticPoint(self.poly,self.x,-self.y,for_factor=self.for_factor)

	def __sub__(self,other):
		other = self.pairwise_check(other)
		return self + (-other)


	def __mul__(self, const:int):
		#like pingala but weird

		mul_st = bin(const)[3:]#drop the '0b' as well as the first bit
		res = self#copy not necessary because all operations create new objects anyway
		for bit in mul_st:
			if bit == '1':
				res = (res + res) + self#r+r+p
			else:
				res = res + res

		return res

	def __hash__(self):
		if self.is_inf:
			return self.p.__hash__()#bad, but it should work
		else:
			return (self.x,self.y).__hash__()


'''
takes parameters as 
y^2 = x^3 + ax + b (mod p)
'''
def EllipticCurve(p,a,b,c=None,**kwargs):
	if c is None:
		c = b
		b = a
		a = 0

	poly = FiniteFieldPoly(p,[1,a,b,c])  #not allowing kwargs here for now
	def get(x=None,y=None):#allow get(), which gives infinity
		return EllipticPoint(poly,x,y,**kwargs)

	return get


'''
factor n using the elliptic method
'''
def elliptic_factor(n, initial_point=None, a=1, verbose=False, ret_mults_to_success=False):
	if initial_point is None:
		initial_point = (ModInteger(0, n), ModInteger(1, n))
	else:
		initial_point = (ModInteger(initial_point[0], n), ModInteger(initial_point[1], n))
	a = ModInteger(a,n)
	#find an EC with that point on it (mod n)
	#this will be
	#y^2 = x^3 + ax + b
	#b = y^2 - (x^3 + ax)
	x0,y0 = initial_point
	b = (y0**2) - ((x0**3) + x0*a)
	E = EllipticCurve(n,a,b,for_factor=True)
	P = E(x0,y0)
	if verbose:
		polystr = str(P.poly)
		print("Polynomial: y² = {}".format(polystr[polystr.index('x'):]))
		print("x₀ = {}, y₀ = {}".format(x0.x,y0.x))
	i = 2
	while not P.factor_canary:
		P *= i
		i += 1
	#we have a winner
	fac0 = ext_eucl_int(P.x.x,n,gcd_only=True)
	fac1 = n//fac0
	if ret_mults_to_success:
		return fac0,fac1,i
	else:
		return fac0,fac1


'''
P also defines the curve itself and the prime 
'''
def ec_el_gamal_enc(m_raw:str,P:EllipticPoint,N:int,target_public_key:EllipticPoint):
	x0s = ascii_to_nums_nomod(m_raw,7)#TODO how do we determine the block size?
	#pad the x0s
	for i in range(len(x0s)):
		x0s[i] *= 100#TODO does this ever change for ascii? (probably not)
	mr = IntegerModRing(P.p,n_is_prime=True)
	#make sure all the f(x)s are square (increment if not) and immediately add to M along with the corresponding y
	Ms = []
	for x0 in x0s:
		fx = P.poly[x0]
		rtfx = fx.sqrt(fast_only=True)  #fast only *should* work
		while rtfx is None:
			#increment x0
			x0 += 1
			if (x0 % 100) == 0:
				raise AttributeError("No square f(x) found for some part of message: {}. Padding size may need to increase.".format(m_raw))
			fx = P.poly[x0]
			rtfx = fx.sqrt(fast_only=True)#fast only *should* work
		#we have a winner
		Ms.append(EllipticPoint(P.poly,x0,rtfx))


	ctxt = []
	for M in Ms:
		#random k
		k = random.randint(1,N)
		S = P*k
		T = M + target_public_key*k
		ctxt.append((S,T))

	return ctxt


def ec_el_gamal_dec(ctxt:List[EllipticPoint],private_key:Union[int,ModInteger]):
	nums = []
	for S,T in ctxt:
		pnum = T - (S*private_key)
		nums.append(pnum.x.x//(100))#divide out the padding

	return nums_to_ascii(nums)


'''
mwahaha

sounds great, doesn't work -- there's way too many complex operations (e.g. EC multiplies, MI inversions) for this to be fast enough
'''
def ec_bday_attack(P:EllipticPoint,B:EllipticPoint,N:int,npoints=None):
	if npoints is None:
		npoints = int(math.sqrt(N) * 1.25)

	#make the lists
	ks = [random.randint(1,N) for _ in range(npoints)]
	ls = [random.randint(1,N) for _ in range(npoints)]
	ls = [l for l in ls if (ext_eucl_int(l,N,gcd_only=True) == 1)]#filter out the ones that won't work anyway (may be a smarter way to do this)
	pks = [P*k for k in ks]
	bls = {B*l:l for l in ls}

	for k,pk in zip(ks,pks):
		if pk in bls:
			#we have a winner
			l = bls[pk]
			#now we have
			#Pk = Bl
			#Pk = bPl
			#Pk = P(bl)
			#so
			#k == bl mod N
			#b == kl^-1 mod N
			b = ModInteger(l,N)**(-1) * k
			return b.x