from crypto_tools.util_upper import *

class ModInteger:
	def __init__(self, x, n, n_is_prime=False, phi_n=None):
		if type(x) == ModInteger:  #if we were created automatically
			self.x = x.x
		else:
			self.x = modi(x, n)
		self.n = n
		if n_is_prime is None:
			self.n_is_prime = isprime(n)
		else:
			self.n_is_prime = n_is_prime

		self.phi_n = phi_n

	def __repr__(self):
		return str(self)

	def __str__(self):
		return "{} (mod {})".format(self.x, self.n)

	def pairwise_check(self, other):
		otherp = other
		if type(other) != ModInteger:
			otherp = ModInteger(other, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)
		elif other.n != self.n:
			raise AttributeError("Pairwise operation on numbers from different moduli")
		return otherp

	def __add__(self, other):
		other = self.pairwise_check(other)
		return ModInteger(self.x+other.x, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)

	def __neg__(self):
		return ModInteger(-self.x, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)

	def __sub__(self, other):
		return self+(-other)

	def __mul__(self, other):
		other = self.pairwise_check(other)
		return ModInteger(self.x*other.x, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)

	def __eq__(self, other):
		other = self.pairwise_check(other)
		return (other.n == self.n) and (other.x == self.x)

	def __ne__(self, other):
		return not self.__eq__(other)

	def __floordiv__(self, other):
		other = self.pairwise_check(other)

		#try running ext eucl to find an inverse, if we find that gcd(other,self) is not 1 then it won't work anyway
		g, (l, _) = ext_eucl_int(other.x,-other.n)
		if g != 1:
			raise ZeroDivisionError(
				"Dividend {} and modulus {} not coprime -- solutions will not be unique".format(other.x, self.n))
		#if we are coprime then l is our inverse

		return self*l

	def __truediv__(self, other):
		return self.__floordiv__(other)

	def __pow__(self, p: int):
		if p == 0:
			return ModInteger(1, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)
		if p < 0:
			return (ModInteger(1, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)/self)**(
				-p)  #a^(-p) = (a^(-1))^p
		return ping(self, p, self.n, n_is_prime=self.n_is_prime, phi_n=self.phi_n)

	def __hash__(self):
		return hash(self.x)

	def __int__(self):
		return self.x

	'''
	give *a* square root modulo n of self, if one exists (if we find one, others may still exist e.g. -a)

	brute-force iff modulo a non prime or p == 1 mod 8
	else with euler-criterion-like methods
	'''
	def sqrt(self, fast_only=False):
		if (self == 0) or (self == 1):  #always squares to itself (this handles the mod 2 case and many others)
			return self
		if fast_only and not self.n_is_prime:  #double-check since this stuff won't work otherwise
			print(
				"WARNING: I don't know if modulus {} is prime and you want me to take fast square roots modulo it.".format(
					self.n))
			self.n_is_prime = isprime(self.n)
		if (self.n%4) == 3:
			sqrt = self**((self.n+1)//4)
			if (sqrt**2) == self:
				return sqrt
		#otherwise sqrt does not exist (-a is square)
		if (self.n%8) == 5:
			sqrt1 = self**((self.n+3)//8)  #option 1
			if (sqrt1**2) == self:
				return sqrt1
			sqrt2 = sqrt1*(ModInteger(2, self.n)**((self.n-1)//4))
			if (sqrt2**2) == self:
				return sqrt2
		#otherwise DNE
		if (not self.n_is_prime) or (fast_only and ((self.n%8) == 1)):
			rt = ModInteger(1, self.n, self.n_is_prime, self.phi_n)
			if self == 0:
				return rt-1  #zero
			while (rt**2 != self) and (rt != 0):
				rt += 1
			if rt == 0:
				return None
			return rt
		return None  #just being explicit (this means root does not exist)


	def multiplicative_order(self):
		#hard to beat brute force
		k = 1
		res = self
		one = ModInteger(1,self.n,self.n_is_prime,self.phi_n)
		while res != one:
			res *= self
			k += 1
		return k

def ext_eucl_int(a,b,gcd_only=False):
	if type(a) == ModInteger:
		a = a.x
	if type(b) == ModInteger:
		b = b.x
	carda = (1,0)
	cardb = (0,1)

	q,r = long_divide(a,b)
	cardc = (carda[0] - (q*cardb[0]),carda[1] - (q*cardb[1]))
	carda = cardb
	cardb = cardc
	a = b
	b = r

	while r != 0:
		q, r = long_divide(a, b)
		cardc = (carda[0]-(q*cardb[0]), carda[1]-(q*cardb[1]))
		carda = cardb
		cardb = cardc
		a = b
		b = r

	if a < 0:
		a = -a
		carda = (-carda[0],-carda[1])

	if gcd_only:
		return a
	else:
		return a,carda


'''
sagelike conversion factory
'''
def IntegerModRing(n,**kwargs):
	def cast_to_mi(x):
		return ModInteger(x,n,**kwargs)
	return cast_to_mi

'''
quickly calculate a^p mod n
'''
def ping(a ,e :int ,n ,**kwargs):
	if ('n_is_prime' not in kwargs) or (kwargs['n_is_prime'] is None):
		kwargs.update({'n_is_prime' :False}  )#this is to deal with infinite recursion problems
	res = ModInteger(a ,n ,**kwargs)
	e_st = bin(e)[3:]
	for bit in e_st:
		if bit == '1':
			#square and multiply by a
			res = res*res*a
		else:
			#just square
			res = res*res

	return res


def primitive_root(p):
	r = IntegerModRing(p)
	m = r(1)
	#factorize phi(p) = p-1
	fpm1 = list(naive_fac(p-1))
	results = [ModInteger(1, p) for _ in range(len(fpm1))]
	while (ModInteger(1, p) in results) and (m.x < p):
		m = m+1
		for i, pf in enumerate(fpm1):
			res_this = m**((p-1)//pf)
			results[i] = res_this

	if ModInteger(1, p) in results:
		raise AttributeError("No primitive root exists modulo {}".format(p))

	return m

def ipret(truth_val,ret_ptrn,primes_to_rootn):
	if ret_ptrn:
		return truth_val,primes_to_rootn
	else:
		return truth_val


'''
miller-rabin test
'''
def mr_is_probprime(n, nbases=None):
	mr = IntegerModRing(n, n_is_prime=False)  #have to specify for the mis that n isn't prime so it doesn't check
	#write n-1 = 2^k * m where m is odd
	#take the binary expansion, k = number of zeroes at the end
	k, m, nbwidth = get_int_2km_fac(n-1, rw=True)
	#choose a base 1 < a < n
	if nbases is None:
		bases = primes(10+(nbwidth>>1))  #good enough for numbers well beyond 64 bits
	else:
		bases = {random.randint(2, n) for _ in range(nbases)}
		while (nbases < n-2) and len(bases) != nbases:
			bases.add(random.randint(2, n))

	for a in bases:
		if a == n:
			return True  #it's prime yo
		#compute b0 = a^m mod n
		b = mr(a)**m
		#if b0 is pm 1 then return probable prime=true

		if (b == 1) or (b == -1):
			continue  #a says probable prime
		#compute b1 = b0^2 mod n
		#if b1 == -1 return pp = true
		#elif b1 == 1 then return prime=false (factor is known)
		#repeat
		m1_seen = False
		for _ in range(k):
			b = b**2
			if b == -1:
				m1_seen = True
				break
			if b == 1:
				return False
		if not m1_seen:
			return False
	#if we get to the end and we didn't get a (-)1 then it's not prime either
	#if we didn't fail yet it's probably prime
	return True

def isprime(n, primes_to_rootn=None, ret_ptrn=False, nbases=None, definite=False):
	mr = IntegerModRing(n, n_is_prime=False)  #obviously n may be prime but this is a code thing
	if n < 2:
		return ipret(False, ret_ptrn, primes_to_rootn)
	if n in [2, 3, 5]:
		return ipret(True, ret_ptrn, primes_to_rootn)
	#do the quick flt check first: a^(n-1) cong 1 mod n if n is prime and a coprime to it
	if not mr_is_probprime(n, nbases=nbases):
		return ipret(False, ret_ptrn, primes_to_rootn)
	elif not definite:
		return ipret(True, ret_ptrn, primes_to_rootn)

	#just do trial division at this point since it might be prime
	rootn = int(np.ceil(np.sqrt(n)))
	if primes_to_rootn is None or rootn > len(primes_to_rootn):
		primes_to_rootn = primes(rootn)
	for p in primes_to_rootn[:rootn]:
		if (n%p) == 0:
			return ipret(False, ret_ptrn, primes_to_rootn)

	return ipret(True, ret_ptrn, primes_to_rootn)


def next_prime(n,check_n=False,**kwargs):
	if check_n:
		if isprime(n,**kwargs):
			return n
	if n < 2:
		return 2
	if n == 2:
		return 3
	p = n + (1 if (n % 2) == 0 else 2)
	p_is_prime = isprime(p,**kwargs)
	while not p_is_prime:
		p += 2
		p_is_prime = isprime(p,**kwargs)
	return p

def prev_prime(n,check_n=False,**kwargs):
	if check_n:
		if isprime(n, **kwargs):
			return n
	if n < 2:
		return None
	if n == 2:
		return 2
	p = n - (1 if (n % 2) == 0 else 2)
	p_is_prime = isprime(p, **kwargs)
	while not p_is_prime:
		p -= 2
		p_is_prime = isprime(p, **kwargs)
	return p


'''
gives a^(-1) mod p where p is prime
'''
def gauss_invert(a,p):
	if a == 1:
		return ModInteger(1,p)
	q,r = long_divide(p,modi(a,p))
	inv_acc = ModInteger(-q,p,n_is_prime=True)
	while r > 1:
		#write p = q*<prev r> + r
		q,r = long_divide(p,r)
		inv_acc = inv_acc * (-q)
	return inv_acc


def naive_fac(n,full=False):
	pl = primes(int(math.sqrt(n))+1)
	f = []
	while not isprime(n):
		#find the lowest prime that divides n and add it to the set, then divide n by that prime
		for p in pl:
			if(n % p == 0):
				f.append(p)
				n //= p
				break
	f.append(n)
	if full:
		return f
	return set(f)


'''
p prime modulus
g primitive root mod p
h in Z/pZ
n is number of trials
'''
def bday_attack_disclog(p, g, h, n=None):
	if n is None:
		n = int(np.sqrt(p)*1.25)
	mr = IntegerModRing(p)
	g = mr(g)
	h = mr(h)

	g_exp_try = list(range(p))
	np.random.shuffle(g_exp_try)
	g_exp_try = g_exp_try[:n]
	gh_exp_try = list(range(p))
	np.random.shuffle(gh_exp_try)
	gh_exp_try = gh_exp_try[:n]

	pg = {(g**exp).x: exp for exp in g_exp_try}
	pgh = {(h*(g**(p-1-exp))).x: exp for exp in gh_exp_try}

	g_exp = None
	gh_exp = None
	for g_res, gh_res in zip(pg, pgh):
		if g_res in pgh:
			g_exp = pg[g_res]
			gh_exp = pgh[g_res]
			break
		if gh_res in pg:
			g_exp = pg[gh_res]
			gh_exp = pgh[gh_res]
			break

	if g_exp is not None:
		return ModInteger(g_exp, p-1)+gh_exp


def ascii_to_nums(s,block_size,phi_n=None):
	mr = IntegerModRing(block_size,phi_n=phi_n)
	i_block_width = 3#technically floor(log10(256))+1 but that's just 3
	s_block_width = (int(math.log10(block_size)) + 1) // i_block_width
	m_nums = []
	for chsq in [s[i:i + s_block_width] for i in range(0, len(s), s_block_width)]:
		nums_this = ''
		for ch in chsq:
			tmp_chord = str(ord(ch))
			tmp_chord = ('0' * (i_block_width - len(tmp_chord))) + tmp_chord
			nums_this += tmp_chord
		m_nums.append(mr(int(nums_this)))
	return m_nums


def el_gamal_enc(p, g, msg, pub_key):
	mr = IntegerModRing(p)
	g = mr(g)

	#convert msg into numbers
	m_nums = ascii_to_nums(msg, p)

	#do all the mults
	rts = []
	for m_num in m_nums:
		#generate a k
		k = np.random.randint(0, p)
		r = g**k
		t = m_num*(mr(pub_key))**(k)

		rts.append((r.x, t.x))

	return rts


'''
m = t*r^(-a)
'''


def el_gamal_dec(p, rts, a, print_chvals=False):
	mr = IntegerModRing(p)
	rts = [(mr(r), mr(t)) for r, t in rts]

	ma = p-1-a
	nums = []
	for r, t in rts:
		m = t*r**(ma)
		#make this into a string, zero-pad, and convert to ascii
		nums.append(m.x)
	return nums_to_ascii(nums, p)