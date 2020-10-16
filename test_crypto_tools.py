from unittest import TestCase
from crypto_tools import *


class Test(TestCase):
	def test_get_qs_factor_base(self):
		#just use the lecture notes example for now
		qsfb, order = get_qs_factor_base(2201, 14, 63)

		assert (1 in qsfb[(5, 1)])
		assert (4 in qsfb[(5, 1)])
		assert (1 in qsfb[(5, 2)])
		assert (24 in qsfb[(5, 2)])
		assert (1 in qsfb[(11, 1)])
		assert (10 in qsfb[(11, 1)])

	# assert (order == [32,16,8,4,2,25,5,11])

	def test_quad_sieve(self):
		facs = quad_sieve(2201, 14, knum=17)
		assert (31 in facs)
		assert (71 in facs)


class TestFiniteFieldPoly(TestCase):

	def test_add(self):
		p = 7
		a = FiniteFieldPoly(p,[4,2,3])
		b = [9,10,11]

		assert((a + b) == [6,5,0])
		b = [1,1]
		assert((a + b) == [4,3,4])

	def test_eq(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = [9, 10, 11]

		assert(a != b)
		b = [7+4,21+2,14+3]
		assert(a == b)

	def test_neg(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = [3, 5, 4]
		assert((-a) == b)

	def test_sub(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = FiniteFieldPoly(p,[9, 10, 11])#2,3,4

		assert((a - b) == [2,-1,-1])
		assert((b-a) == [-2,1,1])

	def test_mul(self):
		p = 7
		a = FiniteFieldPoly(p, [5, 6, 1])
		b = FiniteFieldPoly(p, [1, 2, 0, 5])

		ab = a*b
		assert(np.all(ab.coef == [5,2,6,6,2,5]))

		b = FiniteFieldPoly(p, [1, 1, 1, 1])

		ab = a*b
		assert(np.all(ab.coef == [5,11,12,12,7,1]))

		b = FiniteFieldPoly(p, [2,1])

		ab = b*a
		assert(np.all(ab.coef == [10,17,8,1]))

	def test_div(self):
		p = 7
		a = FiniteFieldPoly(p, [1, 2, 0, 5])
		b = FiniteFieldPoly(p, [5, 6, 1])

		q,r = a//b
		qb = q*b
		qbr = qb + r
		assert(qbr == a)

		p = 2
		a = FiniteFieldPoly(p, [1, -2, 0, -4])
		b = FiniteFieldPoly(p, [1, -3])

		q, r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)
		assert(q == [1,1,3])
		assert(r == [5])

		p = 5
		a = FiniteFieldPoly(p, [3, 0])
		b = FiniteFieldPoly(p, [2])

		q, r = a//b
		qb = q*b
		qbr = qb+r
		assert(qbr == a)
		assert(q == [4,0])
		assert(r == 0)
		
		
	def test_ext_eucl(self):
		p = 5
		a = FiniteFieldPoly(p, [3,0,4,1])
		b = FiniteFieldPoly(p, [2,2,2])
		
		gcd,(s,t) = FFP_ext_eucl(a,b)
		assert(gcd == 2)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert(sm == gcd)

		p = 3
		a = FiniteFieldPoly(p, [1, 0, 0, 1, 1])
		b = FiniteFieldPoly(p, [1, 0, 1])

		gcd, (s, t) = FFP_ext_eucl(a, b)
		assert (gcd == 2)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert (sm == gcd)