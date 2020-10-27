from unittest import TestCase
from crypto_tools import *


class TestModInteger(TestCase):
	def test_sqrt(self):
		assert (ModInteger(0, 2).sqrt(fast_only=True) == 0)
		assert (ModInteger(1, 2).sqrt(fast_only=True) == 1)
		p = 103  #3 mod 4
		for x in range(p):
			x = ModInteger(x, p, n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		while (p%8) != 5:
			p = next_prime(p, check_n=False)
		#5 mod 8
		for x in range(p):
			x = ModInteger(x, p, n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		while (p%8) != 1:
			p = next_prime(p, check_n=False)
		#1 mod 8
		for x in range(p):
			x = ModInteger(x, p, n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		#composite
		n = 80
		for x in range(n):
			x = ModInteger(x, n)
			x2 = x**2
			res = x2.sqrt()
			assert ((res**2) == x2)  #the plus/minus x thing is only true mod a prime


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


class TestFiniteFields(TestCase):

	def test_add(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = [9, 10, 11]

		assert ((a+b) == [6, 5, 0])
		b = [1, 1]
		assert ((a+b) == [4, 3, 4])

	def test_eq(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = [9, 10, 11]

		assert (a != b)
		b = [7+4, 21+2, 14+3]
		assert (a == b)

	def test_neg(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = [3, 5, 4]
		assert ((-a) == b)

	def test_sub(self):
		p = 7
		a = FiniteFieldPoly(p, [4, 2, 3])
		b = FiniteFieldPoly(p, [9, 10, 11])  #2,3,4

		assert ((a-b) == [2, -1, -1])
		assert ((b-a) == [-2, 1, 1])

	def test_mul(self):
		p = 7
		a = FiniteFieldPoly(p, [5, 6, 1])
		b = FiniteFieldPoly(p, [1, 2, 0, 5])

		ab = a*b
		assert (np.all(ab.coef == [5, 2, 6, 6, 2, 5]))

		b = FiniteFieldPoly(p, [1, 1, 1, 1])

		ab = a*b
		assert (np.all(ab.coef == [5, 11, 12, 12, 7, 1]))

		b = FiniteFieldPoly(p, [2, 1])

		ab = b*a
		assert (np.all(ab.coef == [10, 17, 8, 1]))

	def test_div(self):
		p = 7
		a = FiniteFieldPoly(p, [1, 2, 0, 5])
		b = FiniteFieldPoly(p, [5, 6, 1])

		q, r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)

		p = 2
		a = FiniteFieldPoly(p, [1, -2, 0, -4])
		b = FiniteFieldPoly(p, [1, -3])

		q, r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)
		assert (q == [1, 1, 3])
		assert (r == [5])

		p = 5
		a = FiniteFieldPoly(p, [3, 0])
		b = FiniteFieldPoly(p, [2])

		q, r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)
		assert (q == [4, 0])
		assert (r == 0)

	def test_lshift(self):
		p = 3
		a = FiniteFieldPolyModM(p, FiniteFieldPoly(p, [1, 0, 1]), [1, 2, 0, 5])
		e = 5

		r = a<<e
		assert (r == 1)

	def test_ext_eucl(self):
		p = 5
		a = FiniteFieldPoly(p, [3, 0, 4, 1])
		b = FiniteFieldPoly(p, [2, 2, 2])

		gcd, (s, t) = FFP_ext_eucl(a, b)
		assert (gcd == 1)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert (sm == gcd)

		p = 3
		a = FiniteFieldPoly(p, [1, 0, 0, 1, 1])
		b = FiniteFieldPoly(p, [1, 0, 1])

		gcd, (s, t) = FFP_ext_eucl(a, b)
		assert (gcd == 1)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert (sm == gcd)

		p = 17
		a = FiniteFieldPoly(p, [8, 13, 9, 0, 16])
		b = FiniteFieldPoly(p, [1, 0, 0, 0, 0, 0, 12, 14])

		gcd = FFP_ext_eucl(a, b, just_gcd=True)
		assert (gcd == 1)

	def test_eval(self):
		p = 5
		a = FiniteFieldPoly(p, [3, 0, 4, 1])
		b = FiniteFieldPoly(p, [2, 2, 2])
		x = 3

		eva = a[x]
		assert (eva == 4)
		evb = b[x]
		assert (evb == 1)

	def test_ping_ff(self):
		p = 5
		a = FiniteFieldPoly(p, [3, 0, 4, 1])
		e = 5

		ae = a**e
		assert (ae == [243, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2464, 0, 0, 0, 0, 1])
		assert (ae == a*a*a*a*a)

		p = 3
		a = FiniteFieldPolyModM(p, FiniteFieldPoly(p, [1, 0, 1]), [2, 2])  #m=x^2+1, a=2x + 2
		e = 5
		ae = a**e
		assert (ae == [1, 1])
		assert (ae == a*a*a*a*a)

		p = 5
		a = FiniteFieldPoly(p, [3, 0, 4, 1])
		e = 4

		ae = a**e
		assert (ae == a*a*a*a)

	def test_poly_is_reducible(self):
		p = 3
		a = FiniteFieldPoly(p, [1, 0, 1])  #x^2 + 1 is irreducible in F3[x]

		assert (not poly_is_reducible(a))

		a = FiniteFieldPoly(p, [1, 0, 2])  #x^2 + 2 = (x+1)(x+2)

		assert (poly_is_reducible(a))

		p = 17
		a = FiniteFieldPoly(p, [1, 0, 0, 0, 0, 0, 12, 14])  #from sage -- 7th degree (irreducible)

		assert (not poly_is_reducible(a))

		a = FiniteFieldPoly(p, [1, 4, 2, 3])*[5, 5, 3]  #clearly not reducible

		assert (poly_is_reducible(a))

	def test_find_irreducible_poly(self):
		p = 17

		a = find_irreducible_poly(p, 4)
		assert (a == [1, 0, 0, 0, 3])  #regression test


class TestEllipticPoint(TestCase):

	def test_add(self):
		p = 5
		poly = FiniteFieldPoly(p, [1, 0, 2, 4])
		P = EllipticPoint(poly, 0, 2)
		Q = EllipticPoint(poly, 2, 4)
		R = EllipticPoint(poly, 4, 4)

		sm = P+Q
		assert (sm == R)
		R = EllipticPoint(poly, 4, 1)
		sm = P+P
		assert (sm == R)

		P = EllipticPoint(poly, 0, 3)
		Q = EllipticPoint(poly, 4, 1)
		R = EllipticPoint(poly, 0, 2)

		sm = P+Q
		assert (sm == R)

		poly = FiniteFieldPoly(7, [1, 0, 1, 1])
		P = EllipticPoint(poly, 0, 1)
		Q = EllipticPoint(poly, 2, 2)
		R = EllipticPoint(poly, 0, 6)

		sm = P+Q
		assert (sm == R)

		E = EllipticCurve(7, 1, 1, 1)
		P = E(1, 2)
		Q = E(2, 1)
		R = E(4, 1)

		sm = P+Q
		assert (sm == R)

	def test_mul(self):
		p = 19
		E = EllipticCurve(p, 13, 2)
		P = E(10, 7)
		R = E(15, 0)

		pr = P*50
		assert (pr == R)

		p = 101
		E = EllipticCurve(p, 49, 22)
		P = E(53, 81)
		R = E(81, 58)

		pr = P*79
		assert (pr == R)

		E = EllipticCurve(11, 1, 1, 1)
		P = E(2, 2)
		R = E()  #inf

		pr = P*70
		assert (pr == R)

	def test_elliptic_factor(self):
		n = 18923
		f0, f1 = elliptic_factor(n)

		assert ({f0, f1} == {127, 149})

		n = 101*103
		f0, f1 = elliptic_factor(n)

		assert ({f0, f1} == {101, 103})

	def test_ec_el_gamal_enc_and_dec(self):
		E = EllipticCurve(123456789101234567891027, 1, 1, 1)
		P = E(3, 11655832467975276266127)
		N = 61728394550949287614731

		privkey = 666
		pubkey = P*privkey

		#try to encipher a simple message
		msg = 'the quick brown fox Jumped OVER the la7y doogz'
		ctxt = ec_el_gamal_enc(msg, P, N, pubkey)

		#deciphering that should give msg
		dec_msg = ec_el_gamal_dec(ctxt, privkey)
		assert (dec_msg == msg)  #YEET


	def test_ec_bday_attack(self):
		E = EllipticCurve(103,1,1,1)
		P = E(33, 86)
		N = 9
		b = 4
		B = P*b

		rb = ec_bday_attack(P,B,N,npoints=10)
		assert(rb == b)
