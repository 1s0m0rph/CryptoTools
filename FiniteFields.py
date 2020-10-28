from crypto_tools.ModularArith import *

FFP_SS_TRANS = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")  #https://codeigo.com/python/printing-subscript-and-superscript


class FiniteFieldPoly:

	def __init__(self, p: int, coef, print_as_latex=False):
		self.p = p
		self.print_as_latex = print_as_latex
		#build c, being careful with typing
		if type(coef) == FiniteFieldPoly:
			self.coef = coef.coef.copy()
		elif type(coef) in [int, ModInteger]:
			self.coef = np.array([ModInteger(coef, p, n_is_prime=True)])
		else:
			coef_tmp = []
			for c in coef:
				coef_tmp.append(ModInteger(c, p, n_is_prime=True))
			while (len(coef_tmp) > 1) and (coef_tmp[0] == 0):
				del coef_tmp[0]  #no leading zeroes so degree is just the length of coef
			self.coef = np.array(coef_tmp)
		self.dgr = len(self.coef)-1

	def ext_low_degree(self, other):
		low_is_self = self.dgr < other.dgr  #need to return this for non-commutative operations (e.g. subtract, divide)
		ldgc, hdgc = (self.coef, other.coef) if low_is_self else (other.coef, self.coef)
		ldgc_ext = np.array([ModInteger(0, self.p, n_is_prime=True) for _ in range(len(hdgc)-len(ldgc))]+list(ldgc))
		return ldgc_ext, hdgc, low_is_self

	def pairwise_check(self, other):
		if type(other) != FiniteFieldPoly:
			other = FiniteFieldPoly(self.p, other, self.print_as_latex)
			return other, self.ext_low_degree(other)  #all checks will succeed this way, so don't do them
		return other, self.ext_low_degree(other)

	def __repr__(self):
		#make fn look pretty
		if self.print_as_latex:
			fnst = "$"
			for i, c in enumerate(self.coef):
				if (c.x != 0) or (len(self.coef) == 1):
					fnst += ('' if ((c.x == 1) and (i != (len(self.coef)-1))) else "{}".format(c.x))+(
						"x" if (i != ((len(self.coef)-1))) else "")+(
								"^[{}]".format(len(self.coef)-i-1) if i < len(self.coef)-2 else "")+" + "

			fnst = fnst[:-3]  #drop last +

			return fnst.replace('[', '{').replace(']', "}")+"$"
		else:
			fnst = ""
			for i, c in enumerate(self.coef):
				if (c.x != 0) or (len(self.coef) == 1):
					fnst += ('' if ((c.x == 1) and (i != (len(self.coef)-1))) else "{}".format(c.x))+(
						"x" if (i != ((len(self.coef)-1))) else "")+(
								"{}".format(len(self.coef)-i-1).translate(FFP_SS_TRANS) if i < len(
									self.coef)-2 else "")+" + "

			fnst = fnst[:-3]  #drop last +

			return "F{}: {}".format(self.p, fnst)

	def __eq__(self, other):
		if type(other) not in [FiniteFieldPoly, int, list, set, np.array]:
			return False
		other, (ldgc_ext, hdgc, _) = self.pairwise_check(other)
		#as long as that passes all that matters is the degrees
		return np.all(ldgc_ext == hdgc)

	def __add__(self, other):
		other, (ldgc_ext, hdgc, _) = self.pairwise_check(other)
		return FiniteFieldPoly(self.p, ldgc_ext+hdgc, self.print_as_latex)

	def __neg__(self):
		#just negate all the coefficients
		return FiniteFieldPoly(self.p, -self.coef, self.print_as_latex)

	def __sub__(self, other):
		other, (ldgc_ext, hdgc, lis) = self.pairwise_check(other)
		if lis:
			return FiniteFieldPoly(self.p, ldgc_ext-hdgc, self.print_as_latex)
		else:
			return FiniteFieldPoly(self.p, hdgc-ldgc_ext, self.print_as_latex)

	def __lshift__(self, other: int):
		return FiniteFieldPoly(self.p, list(self.coef)+[0 for _ in range(other)], self.print_as_latex)

	def __rshift__(self, other: int):
		return FiniteFieldPoly(self.p, list(self.coef[:-other]), self.print_as_latex)

	'''
	happy happy fun times
	hello karatsuba my old friend

	general idea for karatsuba multiplication:
		want to multiply a*b
		write
			a = a1*B^m + a0
			b = b1*B^m + b0
		for m < min(degree a, degree b)
		then
			ab = (a1B^m + a0)(b1B^m + b0)
			   = c2B^2m + c1B^m + c0
			where
				c2 = a1b1
				c1 = a1b0 + a0b1
				c0 = a0b0
		for all of this use B = x (the parameter of the polynomials)
	'''

	def __mul__(self, other):
		if (type(other) in [ModInteger, int]):
			#special case for scaling
			return FiniteFieldPoly(self.p, self.coef*other, self.print_as_latex)
		other, (ldgc_ext, hdgc, _) = self.pairwise_check(other)
		if len(other.coef) == 1:
			return FiniteFieldPoly(self.p, (self*other.coef[0]).coef, self.print_as_latex)
		if len(self.coef) == 1:
			return FiniteFieldPoly(self.p, (other*self.coef[0]).coef, self.print_as_latex)
		a = hdgc
		b = ldgc_ext
		if len(a) == 1:
			return FiniteFieldPoly(self.p, [a[0]*b[0]], self.print_as_latex)
		#otherwise reduce with karatsuba
		m = len(a)-1
		a0 = FiniteFieldPoly(self.p, a[1:], self.print_as_latex)  #drop the first one since that is a1
		a1 = a[0]
		b0 = FiniteFieldPoly(self.p, b[1:], self.print_as_latex)
		b1 = b[0]
		z2 = FiniteFieldPoly(self.p, a1*b1, self.print_as_latex)  #cheap and easy because a1 and b1 are just constants
		z2scale = z2<<(2*m)  #z2 * B^(2m)
		z0 = a0*b0
		z1 = (a0-a1)*((-b0)+b1)+z2+z0  #this is not infinite since we've reduced the length for a1 and b1
		#ordering for z1 is necessary so we don't have to make more FFPs here
		z1scale = z1<<m
		return z2scale+z1scale+z0

	'''
	poly long divide
	'''

	def __truediv__(self, other):
		#TODO can we use newton-raphson or something similar? (TC: O(n) mults for dividend is degree n, divisor is degree const)
		#this might require generalizing these polys somewhat to have negative exponents
		other, (ldgc_ext, hdgc, lis) = self.pairwise_check(other)
		if lis:
			a = ldgc_ext
			b = hdgc
		else:
			a = hdgc
			b = ldgc_ext

		lcdmi = other.coef[0]**(-1)  #leading coef of the divisor (multiplicative inverse mod p)
		q = FiniteFieldPoly(self.p, 0, self.print_as_latex)
		r = self
		while (r.dgr >= other.dgr) and (r != 0):
			#take the highest-order term (coefficient) from a
			c = r.coef[0]*lcdmi  #MI division (using multiplicative inverse mod p)
			#multiply by 1<<i to obtain the add to q
			aq = (FiniteFieldPoly(self.p, 1, self.print_as_latex)<<(r.dgr-other.dgr))*c
			q += aq
			#multiply by b to obtain this quotient
			qt = other*aq
			#subtract from previous remainder to obtain the next remainder
			r -= qt

		return q, r

	def __reversed__(self):
		return FiniteFieldPoly(self.p, reversed(self.coef), self.print_as_latex)

	'''
	actually evaluate this polynomial

	if the polynomial is called p then evaluating at x is done by
		p[x]
	'''

	def __getitem__(self, x):
		x = ModInteger(x, self.p, n_is_prime=True)
		res = ModInteger(0, self.p, n_is_prime=True)
		for i, c in enumerate(reversed(self.coef)):
			res += c*(x**i)
		return res

	'''
	newton-raphson. update function looks like this:
		xi+1 = xi(2 - a*xi)
	to compute a^(-1) = 1/a

	we're going to do this by using a slight modification on the finite field polynomial idea: allow negative exponents
	that is, we can write x + x^(-1) + 5x^(-3) as <1,0>.<1,0,1>
	for the sake of multiplication when all the exponents are negative, it's okay to treat it as though they were all positive because
		we don't combine unlike terms with addition anyway
		multiplication of (1/x^i)(1/(x^j) is just (1/x^(i+j)) = x^(-i-j)
	when there are positive and negatives involved, we can (at least for the sake of division) just shift both the inverse result and the dividend by enough that the inverse result has only positive exponents (then shift back to get the quotient to the left of the decimal and the remainder to the right)
	thus, this function returns the coefficient list for 1/a, in the opposite order from what they might normally be considered to be
	'''

	def NR_invert(self):
		d = FiniteFieldPoly(self.p, reversed([0]+list(self.coef)),
							self.print_as_latex)  #now we're thinking of this as though its highest term actually has exponent -1 (which acts like 1) technically there's also an x^0 term but since that didn't exist originally we initialize it as zero
		x = FiniteFieldPoly(self.p, 1, self.print_as_latex)
		while x.dgr <= self.dgr:  #this should work as an upper bound since we'll have remainder with degree not larger than ourself
			#do the NR update
			x = x*((-d*x)+2)

		return reversed(x)  #remember to unreverse the coefficients

	def __floordiv__(self, other):
		#TODO can we use newton-raphson or something similar? (TC: O(n) mults for dividend is degree n, divisor is degree const)
		#this might require generalizing these polys somewhat to have negative exponents
		#removed for now -- can't get it to work
		# other, (ldgc_ext, hdgc, lis) = self.pairwise_check(other)
		# if lis:
		# 	a = ldgc_ext
		# 	b = hdgc
		# else:
		# 	a = hdgc
		# 	b = ldgc_ext
		#
		# o_inv = other.NR_invert()
		# #left-shift ourselves to match
		# lss = self << (o_inv.dgr - other.dgr)
		# #multiply
		# resbig = lss * o_inv
		# #last few (reverse order) are remainder
		# r = FiniteFieldPoly(self.p,reversed(resbig.coef[-(self.dgr):-1]))
		# #right-shift for quotient
		# q = resbig >> (2*o_inv.dgr - other.dgr + 1)
		return self/other

	def __pow__(self, power):
		return ping_FF(self, power)


'''
solve the the linear diophantine equation in Fp[x]:
	a(x)s(x) + b(x)t(x) = gcd(a(x),b(x))
'''


def FFP_ext_eucl(a, b, just_gcd=False):
	assert (a.p == b.p)
	ff = FiniteField(a.p)
	acard = (ff(1), ff(0))  #how do we produce a?
	bcard = (ff(0), ff(1))  #how do we produce b?

	q, r = a/b
	ccard = (acard[0]-q*bcard[0], acard[1]-q*bcard[1])
	acard = bcard
	bcard = ccard
	a = b
	b = r

	while r != 0:
		q, r = a/b
		ccard = (acard[0]-q*bcard[0], acard[1]-q*bcard[1])
		acard = bcard
		bcard = ccard
		a = b
		b = r

	if just_gcd:
		#special case for degree zero, but this should probably also exist in general...
		if a.dgr == 0:
			return ext_eucl_int(a.coef[0],a.p,gcd_only=True)
		return a
	else:
		if a.dgr == 0:
			g, (l, _) = ext_eucl_int(a.coef[0].x, a.p)
			if g != 1:
				return g, acard  #acard is technically irrelevant now
			else:
				return ext_eucl_int(a.coef[0], a.p, gcd_only=True), (acard[0]*l, acard[1]*l)  #divide out the a0 term
		return a, acard


def FiniteField(p, **kwargs):
	def gffp(coef):
		return FiniteFieldPoly(p, coef, **kwargs)

	return gffp


'''
Give me natural representatives for a(x) modulo n(x)
'''


def mod_poly(a: FiniteFieldPoly, n: FiniteFieldPoly):
	#TODO: faster algorithm?
	_, r = a/n
	return r


def iterate_poly_coef(poly, extend=False, **kwargs):
	if type(poly) == FiniteFieldPolyModM:
		f = FiniteFieldModM(poly.p, poly.m, **kwargs)
		poly = poly.poly  #want a cracker
	else:
		f = FiniteField(poly.p, **kwargs)
	ncoef = []
	rp = list(reversed(poly.coef))
	for i, c in enumerate(rp):  #big endian isn't strictly necessary
		ncoef.append(c+1)
		if c != -1:
			#add the rest as-is
			for cc in rp[i+1:]:
				ncoef.append(cc)
			break
		elif extend and (i == (len(rp)-1)):
			ncoef.append(ModInteger(1, poly.p))
	return f(reversed(ncoef))


def poly_is_reducible_23(poly):
	#check for roots (these imply reducibility if degree is 2 or 3)
	for r in range(poly.p):
		if poly[r] == 0:
			return True
	return False


'''
rabin irreducibility test
https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
"a polynomial f of degree k in Fq[x] is irreducible iff gcd(f,x^(q^(n_i)) - x) = 1 where n_i is n/p_i where p_i is the ith prime factor of n"
'''


def poly_is_reducible(poly, n_fac_div=None):
	if n_fac_div is None:
		n_fac = naive_fac(poly.dgr)
		n_fac_div = [poly.dgr//pi for pi in n_fac]

	f = FiniteFieldModM(poly.p, poly)

	for fdiv in n_fac_div:
		#calculate the exponent
		e = poly.p**fdiv
		#calculate the base mod f
		h = f(1)<<e  #this is now the polynomial x^(p^(n_i))
		h -= [1, 0]  #h is now the polynomial x^(p^(n_i)) - x
		if ext_eucl_int(h.poly, poly, gcd_only=True) != 1:
			return True
	return False


def find_irreducible_poly(p, min_dgr, **kwargs):
	f = FiniteField(p, **kwargs)
	#"sieve" until we get a hole at at least min_dgr
	if min_dgr == 0:
		return f(0)
	elif min_dgr == 1:
		return f([1, 0])
	elif min_dgr in [2, 3]:
		#first option will not be irreducible so start with [1,0...,1]
		coef_iter = [ModInteger(x, p) for x in [1]+[0]*(min_dgr-1)+[1]]
		p_iter = f(coef_iter)
		while poly_is_reducible_23(p_iter):
			p_iter = iterate_poly_coef(p_iter, extend=False, **kwargs)
		return p_iter

	#so min dgr is at least 4
	#just do trial division -_-
	coef_iter = [ModInteger(x, p) for x in [1]+[0]*(min_dgr-1)+[1]]
	p_iter = f(coef_iter)
	n_fac = naive_fac(p_iter.dgr)
	n_fac_div = [p_iter.dgr//pi for pi in n_fac]
	prev_dgr = p_iter.dgr
	while poly_is_reducible(p_iter, n_fac_div=n_fac_div):
		p_iter = iterate_poly_coef(p_iter, extend=True, **kwargs)
		if p_iter.dgr > prev_dgr:
			prev_dgr = p_iter.dgr  #lol... "p_iter degree"
			n_fac = naive_fac(p_iter.dgr)
			n_fac_div = [p_iter.dgr//pi for pi in n_fac]

	return p_iter


"""
Finite field polynomial modulo m(x) 
"""


class FiniteFieldPolyModM:

	def __init__(self, p: int, m: FiniteFieldPoly, coef, print_as_latex=False):
		assert (m.p == p)
		self.m = m  #we'll take it for granted that this is irreducible for now TODO check & make generator?
		self.p = p
		self.print_as_latex = print_as_latex
		if type(coef) == FiniteFieldPoly:
			self.poly = coef
		else:
			self.poly = FiniteFieldPoly(self.p, coef, self.print_as_latex)
		#now mod that with m
		self.poly = mod_poly(self.poly, self.m)
		self.inv = None  #cache this

	def __repr__(self):
		if self.print_as_latex:
			return "{}".format(self.poly)
		else:
			m_fmt = "{}".format(self.m)
			m_fmt = m_fmt[m_fmt.index(' '):]  #drop the "F<p>: " part from m only
			return "{}  (mod {})".format(self.poly, m_fmt)

	def pairwise_check(self, other):
		if type(other) != FiniteFieldPolyModM:
			other = FiniteFieldPolyModM(self.p, self.m, other, self.print_as_latex)
			return other
		if other.p != self.p:
			raise AttributeError("Polynomials are in different fields")
		if other.m != self.m:
			raise AttributeError("Polynomials are modulo different moduli")
		return other

	def bin_op_std(self, other, op):
		other = self.pairwise_check(other)
		return FiniteFieldPolyModM(self.p, self.m, op(self.poly, other.poly), self.print_as_latex)

	def __eq__(self, other):
		other = self.pairwise_check(other)
		return self.poly == other.poly

	def __add__(self, other):
		return self.bin_op_std(other, lambda x, y: x+y)

	def __sub__(self, other):
		return self.bin_op_std(other, lambda x, y: x-y)

	def __neg__(self):
		return FiniteFieldPolyModM(self.p, self.m, -self.poly, self.print_as_latex)

	def __mul__(self, other):
		return self.bin_op_std(other, lambda x, y: x*y)

	def __pow__(self, power):
		return ping_FF(self, power)

	def __lshift__(self, other: int):
		#keep going until the degree gets bigger than m, then mod, mod again at the end
		res = FiniteFieldPoly(self.p, self.poly, self.print_as_latex)  #copy
		k = 0
		while k < other:
			shift_amt = min(other-k, self.m.dgr-res.dgr)  #this is the exact required amount
			res <<= shift_amt
			res = mod_poly(res, self.m)
			k += shift_amt
		return FiniteFieldPolyModM(self.p, self.m, res, self.print_as_latex)

	'''
	find the multiplicative inverse of self mod m (ext eucl)
	equivalent to solving the linear diophantine eq:
		a(x)*a^-1(x) - m(x)k(x) = 1
	'''

	def get_inv(self):
		if self.inv is not None:
			return self.inv
		g, (self_inv, _) = FFP_ext_eucl(self.poly, -self.m)
		if g != 1:
			raise AttributeError("GCD of modulus and polynomial is not <1>, no inverse exists")
		si = FiniteFieldPolyModM(self.p, self.m, self_inv, self.print_as_latex)
		si.inv = self
		self.inv = si
		return si

	def __floordiv__(self, other):
		other = self.pairwise_check(other)
		return self*other.get_inv()

	def __truediv__(self, other):
		return self.__floordiv__(other)


def FiniteFieldModM(p, m, **kwargs):
	def get(coef):
		return FiniteFieldPolyModM(p, m, coef, **kwargs)

	return get


'''
pingala implementation for finite fields
'''


def ping_FF(fp: Union[FiniteFieldPoly, FiniteFieldPolyModM], e: int):
	if type(fp) == FiniteFieldPoly:
		poly = fp
		res = FiniteFieldPoly(poly.p, poly.coef, poly.print_as_latex)
	else:
		poly = fp.poly
		res = FiniteFieldPolyModM(fp.p, fp.m, fp.poly, fp.print_as_latex)

	eb_st = bin(e)[3:]  #drop the '0b' as well as the first bit since it's always a 1
	for bit in eb_st:
		if bit == '1':
			res = res*res*poly
		else:
			res = res*res

	return res


'''
Finite Field Polynomial Mod M multiplication table in latex
'''


def FFPMM_mult_table_tex(p, m):
	#first find the elements of F_p[x]/(m(x))
	residues = [FiniteFieldPolyModM(p, m, 0, print_as_latex=True)]
	#build using polynomial iteration from zero
	for _ in range((p**(m.dgr))-1):  #already did zero
		residues.append(iterate_poly_coef(residues[-1], extend=True, print_as_latex=True))

	#now build the actual table
	table = [[None for _ in range(len(residues))] for _ in range(len(residues))]
	for i, r1 in enumerate(residues):
		for j in range(i, len(residues)):
			r2 = residues[j]
			table[i][j] = table[j][i] = r1*r2

	#print (latexified)
	print(r"\begin{center}")
	print(r"\begin{tabular}{c||"+'c|'*len(residues)+"}")
	print(r"$\times$ & ", end='')
	for item in residues[:-1]:
		print("{} & ".format(item), end='')
	print("{}\\\\\\hline\\hline".format(residues[-1]))

	for i, row in enumerate(table):
		print("{} & ".format(residues[i]), end='')
		for item in row[:-1]:
			print("{} & ".format(item), end='')
		print("{}\\\\\\hline".format(row[-1]))
	print(r"\end{tabular}")
	print(r"\end{center}")