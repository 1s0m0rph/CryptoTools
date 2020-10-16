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
		ax = FiniteFieldPoly(p,[4,2,3])
		bx = [9,10,11]

		assert((ax + bx) == [6,5,0])
		bx = [1,1]
		assert((ax + bx) == [4,3,4])

	def test_eq(self):
		p = 7
		ax = FiniteFieldPoly(p, [4, 2, 3])
		bx = [9, 10, 11]

		assert(ax != bx)
		bx = [7+4,21+2,14+3]
		assert(ax == bx)

	def test_neg(self):
		p = 7
		ax = FiniteFieldPoly(p, [4, 2, 3])
		bx = [3, 5, 4]
		assert((-ax) == bx)

	def test_sub(self):
		p = 7
		ax = FiniteFieldPoly(p, [4, 2, 3])
		bx = FiniteFieldPoly(p,[9, 10, 11])#2,3,4

		assert((ax - bx) == [2,-1,-1])
		assert((bx-ax) == [-2,1,1])

	def test_mul(self):
		p = 7
		ax = FiniteFieldPoly(p, [5, 6, 1])
		bx = FiniteFieldPoly(p, [1, 2, 0, 5])

		axbx = ax*bx
		assert(np.all(axbx.coef == [5,2,6,6,2,5]))

		bx = FiniteFieldPoly(p, [1, 1, 1, 1])

		axbx = ax*bx
		assert(np.all(axbx.coef == [5,11,12,12,7,1]))

		bx = FiniteFieldPoly(p, [2,1])

		axbx = bx*ax
		assert(np.all(axbx.coef == [10,17,8,1]))

	def test_div(self):
		p = 7
		ax = FiniteFieldPoly(p, [1, 2, 0, 5])
		bx = FiniteFieldPoly(p, [5, 6, 1])

		q,r = ax//bx
		qbx = q*bx
		qbxr = qbx + r
		assert(qbxr == ax)

		p = 2
		ax = FiniteFieldPoly(p, [1, -2, 0, -4])
		bx = FiniteFieldPoly(p, [1, -3])

		q, r = ax//bx
		qbx = q*bx
		qbxr = qbx+r
		assert (qbxr == ax)
		assert(q == [1,1,3])
		assert(r == [5])