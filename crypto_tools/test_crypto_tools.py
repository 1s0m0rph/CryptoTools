from unittest import TestCase
from crypto_tools.interface import *


class TestModInteger(TestCase):
	def test_sqrt(self):
		assert (ModInteger(0,2).sqrt(fast_only=True) == 0)
		assert (ModInteger(1,2).sqrt(fast_only=True) == 1)
		p = 103  #3 mod 4
		for x in range(p):
			x = ModInteger(x,p,n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		while (p%8) != 5:
			p = next_prime(p,check_n=False)
		#5 mod 8
		for x in range(p):
			x = ModInteger(x,p,n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		while (p%8) != 1:
			p = next_prime(p,check_n=False)
		#1 mod 8
		for x in range(p):
			x = ModInteger(x,p,n_is_prime=True)
			res = (x**2).sqrt(fast_only=True)
			assert ((res == x) or (res == -x))

		#composite
		n = 80
		for x in range(n):
			x = ModInteger(x,n)
			x2 = x**2
			res = x2.sqrt()
			assert ((res**2) == x2)  #the plus/minus x thing is only true mod a prime


	def test_mult_order(self):
		n = 7
		R = IntegerModRing(n)
		a = R(3)

		assert(a.multiplicative_order() == 6)

		a = R(2)

		assert(a.multiplicative_order() == 3)

		n = 1001
		R = IntegerModRing(n)
		a = R(101)

		assert(a.multiplicative_order() == 30)



class TestQuadSieve(TestCase):
	def test_get_qs_factor_base(self):
		#just use the lecture notes example for now
		qsfb,order = get_qs_factor_base(2201,14,63)

		assert (1 in qsfb[(5,1)])
		assert (4 in qsfb[(5,1)])
		assert (1 in qsfb[(5,2)])
		assert (24 in qsfb[(5,2)])
		assert (1 in qsfb[(11,1)])
		assert (10 in qsfb[(11,1)])

	# assert (order == [32,16,8,4,2,25,5,11])

	def test_quad_sieve(self):
		facs = quad_sieve(2201,14,knum=17)
		assert (31 in facs)
		assert (71 in facs)


class TestFiniteFields(TestCase):

	def test_add(self):
		p = 7
		a = FiniteFieldPoly(p,[4,2,3])
		b = [9,10,11]

		assert ((a+b) == [6,5,0])
		b = [1,1]
		assert ((a+b) == [4,3,4])

	def test_eq(self):
		p = 7
		a = FiniteFieldPoly(p,[4,2,3])
		b = [9,10,11]

		assert (a != b)
		b = [7+4,21+2,14+3]
		assert (a == b)

	def test_neg(self):
		p = 7
		a = FiniteFieldPoly(p,[4,2,3])
		b = [3,5,4]
		assert ((-a) == b)

	def test_sub(self):
		p = 7
		a = FiniteFieldPoly(p,[4,2,3])
		b = FiniteFieldPoly(p,[9,10,11])  #2,3,4

		assert ((a-b) == [2,-1,-1])
		assert ((b-a) == [-2,1,1])

	def test_mul(self):
		p = 7
		a = FiniteFieldPoly(p,[5,6,1])
		b = FiniteFieldPoly(p,[1,2,0,5])

		ab = a*b
		assert (np.all(ab.coef == [5,2,6,6,2,5]))

		b = FiniteFieldPoly(p,[1,1,1,1])

		ab = a*b
		assert (np.all(ab.coef == [5,11,12,12,7,1]))

		b = FiniteFieldPoly(p,[2,1])

		ab = b*a
		assert (np.all(ab.coef == [10,17,8,1]))


		#test case that cropped up in lattices
		p = 101
		m = FiniteFieldPoly(p,[1,0,0,0,1])
		F = FiniteFieldModM(p,m)
		a = F([83,23,51,77])
		b = F([96,97,26,74])
		r = F([100,100,0])
		e1 = F([1,0,0])
		msg = 3
		k = 20
		e2 = F([100,1,0])

		v_exp = F([27,75,6,5])#SHOULD be this (from notes)
		v = (a*r) + e1
		assert(v == v_exp)

		w_exp = F([79,0,23,51])
		w = (b*r) + e2 + (msg*k)
		assert(w_exp == w)

	def test_div(self):
		p = 7
		a = FiniteFieldPoly(p,[1,2,0,5])
		b = FiniteFieldPoly(p,[5,6,1])

		q,r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)

		p = 2
		a = FiniteFieldPoly(p,[1,-2,0,-4])
		b = FiniteFieldPoly(p,[1,-3])

		q,r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)
		assert (q == [1,1,3])
		assert (r == [5])

		p = 5
		a = FiniteFieldPoly(p,[3,0])
		b = FiniteFieldPoly(p,[2])

		q,r = a//b
		qb = q*b
		qbr = qb+r
		assert (qbr == a)
		assert (q == [4,0])
		assert (r == 0)

	def test_lshift(self):
		p = 3
		a = FiniteFieldPolyModM(p,FiniteFieldPoly(p,[1,0,1]),[1,2,0,5])
		e = 5

		r = a<<e
		assert (r == 1)

	def test_ext_eucl(self):
		p = 5
		a = FiniteFieldPoly(p,[3,0,4,1])
		b = FiniteFieldPoly(p,[2,2,2])

		gcd,(s,t) = FFP_ext_eucl(a,b)
		assert (gcd == 1)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert (sm == gcd)

		p = 3
		a = FiniteFieldPoly(p,[1,0,0,1,1])
		b = FiniteFieldPoly(p,[1,0,1])

		gcd,(s,t) = FFP_ext_eucl(a,b)
		assert (gcd == 1)
		ats = a*s
		btt = b*t
		sm = ats+btt
		assert (sm == gcd)

		p = 17
		a = FiniteFieldPoly(p,[8,13,9,0,16])
		b = FiniteFieldPoly(p,[1,0,0,0,0,0,12,14])

		gcd = FFP_ext_eucl(a,b,gcd_only=True)
		assert (gcd == 1)

	def test_eval(self):
		p = 5
		a = FiniteFieldPoly(p,[3,0,4,1])
		b = FiniteFieldPoly(p,[2,2,2])
		x = 3

		eva = a[x]
		assert (eva == 4)
		evb = b[x]
		assert (evb == 1)

	def test_ping_ff(self):
		p = 5
		a = FiniteFieldPoly(p,[3,0,4,1])
		e = 5

		ae = a**e
		assert (ae == [243,0,0,0,0,0,0,0,0,0,2464,0,0,0,0,1])
		assert (ae == a*a*a*a*a)

		p = 3
		a = FiniteFieldPolyModM(p,FiniteFieldPoly(p,[1,0,1]),[2,2])  #m=x^2+1, a=2x + 2
		e = 5
		ae = a**e
		assert (ae == [1,1])
		assert (ae == a*a*a*a*a)

		p = 5
		a = FiniteFieldPoly(p,[3,0,4,1])
		e = 4

		ae = a**e
		assert (ae == a*a*a*a)

	def test_poly_is_reducible(self):
		p = 3
		a = FiniteFieldPoly(p,[1,0,1])  #x^2 + 1 is irreducible in F3[x]

		assert (not poly_is_reducible(a))

		a = FiniteFieldPoly(p,[1,0,2])  #x^2 + 2 = (x+1)(x+2)

		assert (poly_is_reducible(a))

		p = 17
		a = FiniteFieldPoly(p,[1,0,0,0,0,0,12,14])  #from sage -- 7th degree (irreducible)

		assert (not poly_is_reducible(a))

		a = FiniteFieldPoly(p,[1,4,2,3])*[5,5,3]  #clearly not reducible

		assert (poly_is_reducible(a))

	def test_find_irreducible_poly(self):
		p = 17

		a = find_irreducible_poly(p,4)
		assert (a == [1,0,0,0,3])  #regression test


class TestEllipticPoint(TestCase):

	def test_add(self):
		p = 5
		poly = FiniteFieldPoly(p,[1,0,2,4])
		P = EllipticPoint(poly,0,2)
		Q = EllipticPoint(poly,2,4)
		R = EllipticPoint(poly,4,4)

		sm = P+Q
		assert (sm == R)
		R = EllipticPoint(poly,4,1)
		sm = P+P
		assert (sm == R)

		P = EllipticPoint(poly,0,3)
		Q = EllipticPoint(poly,4,1)
		R = EllipticPoint(poly,0,2)

		sm = P+Q
		assert (sm == R)

		poly = FiniteFieldPoly(7,[1,0,1,1])
		P = EllipticPoint(poly,0,1)
		Q = EllipticPoint(poly,2,2)
		R = EllipticPoint(poly,0,6)

		sm = P+Q
		assert (sm == R)

		E = EllipticCurve(7,1,1,1)
		P = E(1,2)
		Q = E(2,1)
		R = E(4,1)

		sm = P+Q
		assert (sm == R)

	def test_mul(self):
		p = 19
		E = EllipticCurve(p,13,2)
		P = E(10,7)
		R = E(15,0)

		pr = P*50
		assert (pr == R)

		p = 101
		E = EllipticCurve(p,49,22)
		P = E(53,81)
		R = E(81,58)

		pr = P*79
		assert (pr == R)

		E = EllipticCurve(11,1,1,1)
		P = E(2,2)
		R = E()  #inf

		pr = P*70
		assert (pr == R)

	def test_elliptic_factor(self):
		n = 18923
		f0,f1 = elliptic_factor(n)

		assert ({f0,f1} == {127,149})

		n = 101*103
		f0,f1 = elliptic_factor(n)

		assert ({f0,f1} == {101,103})

	def test_ec_el_gamal_enc_and_dec(self):
		E = EllipticCurve(123456789101234567891027,1,1,1)
		P = E(3,11655832467975276266127)
		N = 61728394550949287614731

		privkey = 666
		pubkey = P*privkey

		#try to encipher a simple message
		msg = 'the quick brown fox Jumped OVER the la7y doogz'
		ctxt = ec_el_gamal_enc(msg,P,N,pubkey)

		#deciphering that should give msg
		dec_msg = ec_el_gamal_dec(ctxt,privkey)
		assert (dec_msg == msg)  #YEET

	def test_ec_bday_attack(self):
		E = EllipticCurve(103,1,1,1)
		P = E(33,86)
		N = 9
		b = 4
		B = P*b

		rb = ec_bday_attack(P,B,N,npoints=10)
		assert (rb == b)

class TestRLWE(TestCase):

	def test_generate_keypair(self):
		#just blind test for now
		#baby example from the notes 2020-11-06
		p = 101
		n = 4
		k = 20
		R = 1
		S = RLWE(p,n,k,R)

		#start by verifying the private key from the actual example
		private_key = FiniteFieldPoly(p,[1,0,0,-1])
		_,a,b =  S.generate_keypair(private_key=private_key)#since this is really random there's not much for expectations

		#do a random one as well and check smallness
		private_key,a,b = S.generate_keypair()

		for c in private_key.poly.coef:
			assert((c.x <= R) or ((-c).x <= R))

		#input stuff
		S = RLWE(p,n,k,R,print_for_pyinput=True)
		s,a,b = S.generate_keypair()
		pass#this is really just for debug mode


	def test_single_message_en_decrypt(self):
		#baby example from the notes 2020-11-06
		p = 101
		n = 4
		k = 20
		R = 1
		S = RLWE(p,n,k,R)

		#keypair from lecture
		s = FiniteFieldPolyModM(p,S.m,[1,0,0,-1])
		a = FiniteFieldPolyModM(p,S.m,[83,23,51,77])
		b = FiniteFieldPolyModM(p,S.m,[96,97,26,74])

		m = 3
		random.seed(1)

		v,w = S.encrypt_single_message_to_key(a,b,m)#this depends on an ephemeral key so it won't be predictable

		m_rec = S.decrypt_single_message(v,w,s)

		assert(m_rec == m)

		#bigger example
		p = 100000000003
		n = 8
		k = 100
		R = 1
		S = RLWE(p,n,k,R)

		#random keypair
		s,a,b =  S.generate_keypair()

		m = 123146852#just anything with less than 10 digits should  work

		v,w = S.encrypt_single_message_to_key(a,b,m)  #this depends on an ephemeral key so it won't be predictable

		m_rec = S.decrypt_single_message(v,w,s)

		assert(m_rec == m)


	def test_str_en_decrypt(self):
		#bigger example
		p = 100000000003
		n = 8
		k = 100
		R = 1
		S = RLWE(p,n,k,R)

		#random keypair
		s,a,b = S.generate_keypair()

		#encrypt a long message to a,b
		msg_orig = "this is the message"
		ctxt = S.encrypt(a,b,msg_orig)

		#decrypt that
		msg_rec = S.decrypt(ctxt,s)

		assert(msg_rec == msg_orig)

		#since this didn't work in practice for some reason...
		F = FiniteFieldModM(p,S.m)
		s = F([-1,-1,1,0,0,0,-1,-1])
		v = F([58354816316,4597505578,27080648777,40612538177,17748688086,76669112944,13884701424,83476093754])
		w = F([16355995873,57085269629,78735228696,14962111776,68677880351,40726763122,56444467681,58393960943])
		msg_orig = 'key'

		ctxt = [(v,w)]
		msg_rec = S.decrypt(ctxt,s)

		assert(msg_rec == msg_orig)

class TestRational(TestCase):
	def test_add(self):
		x = Rational(11,8)
		y = Rational(53,200)

		s = x + y
		assert(s == Rational(41,25))

	def test_sub(self):
		x = Rational(11,8)
		y = Rational(53,200)

		d = x-y
		assert (d == Rational(111,100))

	def test_inv(self):
		x = Rational(11,8)

		ix = ~x
		assert(ix == Rational(8,11))

	def test_div(self):
		x = Rational(11,8)
		y = Rational(53,200)

		d = x/y
		assert (d == Rational(275,53))


	def test_continued_frac_convergents(self):
		x = Rational(17,39)

		convs = continued_frac_convergents(x)
		econvs = [Rational(0,1),
				  Rational(1,2),
				  Rational(3,7),
				  Rational(7,16),
				  Rational(17,39)]

		assert(convs == econvs)

		x = Rational(275,53)

		convs = continued_frac_convergents(x)
		econvs = [Rational(5,1),
				  Rational(26,5),
				  Rational(83,16),
				  Rational(275,53)]

		assert (convs == econvs)

	def test_rp(self):
		x = rp('1/2')

		assert(x.a == 1)
		assert(x.b == 2)

		x = rp('12')

		assert(x.a == 12)
		assert(x.b == 1)

		x = rp('2/2')

		assert(x.a == 1)
		assert(x.b == 1)

	def test_p_expansion(self):
		x = Rational(1,3)
		x_i_exp = 0
		x_ni_exp = [3,3,3,3,3,3,3,3,3,3]

		x_i,x_ni = x.p_expansion(10)
		assert (x_i_exp == x_i)
		assert (x_ni_exp == x_ni)

		x = Rational(15,64)
		x_i_exp = 0
		x_ni_exp = [2,3,4,3,7,5,0,0,0,0]

		x_i,x_ni = x.p_expansion(10)
		assert (x_i_exp == x_i)
		assert (x_ni_exp == x_ni)

		x = Rational(303,49)
		x_i_exp = 6
		x_ni_exp = [1,8,3,6,7,3,4,6,9,3,8,7,7]

		x_i,x_ni = x.p_expansion(13)
		assert (x_i_exp == x_i)
		assert (x_ni_exp == x_ni)