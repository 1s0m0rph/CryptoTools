"""
BASIC NUMBER THEORY CODE
"""

import numpy as np
import math
import random

PRIMES_1K = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997]

ENGL_FREQ_DICT = {'A':0.082,
				'B':0.015,
				'C':0.028,
				'D':0.043,
				'E':0.13,
				'F':0.022,
				'G':0.02,
				'H':0.061,
				'I':0.07,
				'J':0.0015,
				'K':0.0077,
				'L':0.04,
				'M':0.024,
				'N':0.067,
				'O':0.075,
				'P':0.019,
				'Q':0.00095,
				'R':0.06,
				'S':0.063,
				'T':0.091,
				'U':0.028,
				'V':0.0098,
				'W':0.024,
				'X':0.0015,
				'Y':0.02,
				'Z':0.00074
}

ENGL_FREQ = [ENGL_FREQ_DICT[k] for k in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ']

def mod(a,n):
	if a < 0:
		return ModInteger(n*(-(a // n)) + a,n)
	return ModInteger(a % n,n)

def modi(a,n):
	if a < 0:
		return n*(-(a // n)) + a
	return a % n

def gcd(a,b,returnSteps=False,modOperation = modi):
	a = abs(a)
	b = abs(b)
	currentA = max(a,b)
	currentB = min(a,b)
	if currentB == 0:
		return currentA
	if(returnSteps):
		steps = [(currentA,currentB)]
	rem = modOperation(currentA,currentB)
	while(rem != 0):
		currentA = abs(currentB)
		currentB = abs(rem)
		rem = modOperation(currentA,currentB)
		if(returnSteps):
			steps.append((currentA,currentB))
	if(returnSteps):
		steps.append((currentB,rem))
		return currentB,steps
	return currentB

def ext_eucl(a,b,returnSteps=False):
	modOperation = modi
	a = abs(a)
	b = abs(b)
	currentA = a
	currentB = b
	cardca = (1,0)
	cardcb = (0,1)
	if currentB == 0:
		return currentA,cardca
	if currentA == 0:
		return currentB,cardcb
	if(returnSteps):
		steps = []
	rem = modOperation(currentA,currentB)
	while(rem != 0):
		k = currentA//currentB
		cardcc = (cardca[0]-k*cardcb[0],cardca[1]-k*cardcb[1])
		if(returnSteps):
			steps.append((currentA,cardca,currentB,cardcb,rem,cardcc,k))
		cardca = cardcb
		cardcb = cardcc
		currentA = abs(currentB)
		currentB = abs(rem)
		rem = modOperation(currentA,currentB)
	if(returnSteps):
		k = int(currentA/currentB)
		cardcc = (cardca[0]-k*cardcb[0],cardca[1]-k*cardcb[1])
		steps.append((currentA,cardca,currentB,cardcb,rem,cardcc,k))
		return currentB,cardcb,steps
	return currentB,cardcb

def primes(n):
	if n <= 1000:
		i = 0
		while (i < len(PRIMES_1K)) and (PRIMES_1K[i] <= n):
			i += 1
		return PRIMES_1K[:i]

	#seive up to n
	r = list(range(2,n+1))
	pi = 0
	while(pi <= len(r)-1):
		jmp = r[pi]
		for j in range(pi+jmp,len(r),jmp):
			r[j] = 0
		pi += 1
		while(pi <= len(r) - 1) and (r[pi] == 0):
			pi += 1
	i = 0
	while(i < len(r)):
		if r[i] == 0:
			del r[i]
		else:
			i += 1
	return r

'''
composite numbers up to n
'''
def composites(n):
	if n < 4:
		return {}
	all_to_n = set(range(4,n))
	p_to_n = set(primes(n))
	return all_to_n - p_to_n

def naive_fac(n):
	pl = primes(int(math.sqrt(n))+1)
	ps = set(pl)
	f = []
	while not isprime(n):
		#find the lowest prime that divides n and add it to the set, then divide n by that prime
		for p in pl:
			if(n % p == 0):
				f.append(p)
				n //= p
				break
	f.append(n)
	return set(f)

def ipret(truth_val,ret_ptrn,primes_to_rootn):
	if ret_ptrn:
		return truth_val,primes_to_rootn
	else:
		return truth_val

'''
miller-rabin test
'''
def mr_is_probprime(n,nbases=None):
	mr = IntegerModRing(n,n_is_prime=False)#have to specify for the mis that n isn't prime so it doesn't check
	#write n-1 = 2^k * m where m is odd
	#take the binary expansion, k = number of zeroes at the end
	nbin = list(reversed(bin(n-1)[2:]))
	nbwidth = len(nbin)
	k = 0
	while (len(nbin) > 0) and (nbin[0] == '0'):
		k += 1
		del nbin[0]
	#remaining part is m
	nbin = ''.join(reversed(nbin))
	m = int(nbin,base=2)
	#choose a base 1 < a < n
	if nbases is None:
		bases = primes(10 + (nbwidth>>1))#good enough for numbers well beyond 64 bits
	else:
		bases = {random.randint(2,n) for _ in range(nbases)}
		while (nbases < n-2) and len(bases) != nbases:
			bases.add(random.randint(2,n))
	
	for a in bases:
		#compute b0 = a^m mod n
		b = mr(a)**m
		#if b0 is pm 1 then return probable prime=true
		
		if (b == 1) or (b == -1):
			continue#a says probable prime
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
	
def isprime(n,primes_to_rootn=None,ret_ptrn=False,nbases=None,definite=False):
	mr = IntegerModRing(n,n_is_prime=False)#obviously n may be prime but this is a code thing
	if n < 2:
		return ipret(False,ret_ptrn,primes_to_rootn)
	if n in [2,3,5]:
		return ipret(True,ret_ptrn,primes_to_rootn)
	#do the quick flt check first: a^(n-1) cong 1 mod n if n is prime and a coprime to it
	if not mr_is_probprime(n,nbases=nbases):
		return ipret(False,ret_ptrn,primes_to_rootn)
	elif not definite:
		return ipret(True,ret_ptrn,primes_to_rootn)
		
	#just do trial division at this point since it might be prime
	rootn = int(np.ceil(np.sqrt(n)))
	if primes_to_rootn is None or rootn > len(primes_to_rootn):
		primes_to_rootn = primes(rootn)
	for p in primes_to_rootn[:rootn]:
		if (n % p) == 0:
			return ipret(False,ret_ptrn,primes_to_rootn)
	
	return ipret(True,ret_ptrn,primes_to_rootn)

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

def primitive_root(p):
	r = IntegerModRing(p)
	m = r(1)
	#factorize phi(p) = p-1
	fpm1 = list(naive_fac(p-1))
	results = [ModInteger(1,p) for _ in range(len(fpm1))]
	while (ModInteger(1,p) in results) and (m.x < p):
		m = m + 1
		for i,pf in enumerate(fpm1):
			res_this = m**((p-1)//pf)
			results[i] = res_this
	
	if ModInteger(1,p) in results:
		raise AttributeError("No primitive root exists modulo {}".format(p))
	
	return m
	

'''
get most significant bit of x
'''
def get_msb(x):
	return len(bin(x))-2

'''
quickly calculate a^p mod n
'''
def ping(a,e:int,n):
	res = ModInteger(a,n,n_is_prime=False)#n may be prime, but we don't want to check
	e_st = bin(e)[3:]
	for bit in e_st:
		if bit == '1':
			#square and multiply by a
			res = res*res*a
		else:
			#just square
			res = res*res
	
	return res

def long_divide(a,b):
	q = a//b
	if (a - q*b) < (b/2):
		r = a - q*b
	else:
		q += 1
		r = q*b - a
		q = -q
	return q,r

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

class ModInteger:
	def __init__(self, x, n, n_is_prime=False, phi_n=None):
		if type(x) == ModInteger:#if we were created automatically
			self.x = x.x
		else:
			self.x = modi(x,n)
		self.n = n
		if n_is_prime is None:
			self.n_is_prime = isprime(n)
		else:
			self.n_is_prime = n_is_prime
		
		self.phi_n = phi_n
		
		
	def __repr__(self):
		return str(self)
	
	def __str__(self):
		return "{} (mod {})".format(self.x,self.n)
	
	def pairwise_check(self,other):
		otherp = other
		if type(other) != ModInteger:
			otherp = ModInteger(other,self.n,self.n_is_prime,self.phi_n)
		elif other.n != self.n:
			raise AttributeError("Pairwise operation on numbers from different moduli")
		return otherp
	
	def __add__(self,other):
		other = self.pairwise_check(other)
		return ModInteger(self.x + other.x,self.n,self.n_is_prime,self.phi_n)
	
	def __neg__(self):
		return ModInteger(-self.x,self.n,self.n_is_prime,self.phi_n)
	
	def __sub__(self,other):
		return self + (-other)
	
	def __mul__(self,other):
		other = self.pairwise_check(other)
		return ModInteger(self.x * other.x,self.n,self.n_is_prime,self.phi_n)
	
	def __eq__(self,other):
		other = self.pairwise_check(other)
		return (other.n == self.n) and (other.x == self.x)
	
	def __ne__(self,other):
		return not self.__eq__(other)
	
	def __floordiv__(self,other):
		other = self.pairwise_check(other)

		#try running ext eucl to find an inverse, if we find that gcd(other,self) is not 1 then it won't work anyway
		g,(l,_) = ext_eucl(other.x,other.n)
		if g != 1:
			raise AttributeError("Dividend and modulus not coprime -- solutions will not be unique")
		#if we are coprime then l is our inverse

		return self * l
		
	def __truediv__(self,other):
		return self.__floordiv__(other)
	
	def __pow__(self,p:int):
		if p == 0:
			return ModInteger(1,self.n,self.n_is_prime,self.phi_n)
		if p < 0:
			return (ModInteger(1,self.n,self.n_is_prime,self.phi_n)/self)**(-p)#a^(-p) = (a^(-1))^p
		return ping(self,p,self.n)

'''
sagelike conversion factory
'''
def IntegerModRing(n,**kwargs):
	def cast_to_mi(x):
		return ModInteger(x,n,**kwargs)
	return cast_to_mi

def b_to_dec(n,bsym):
	#invert the bsymbols
	bsym_inv = {ch:i for i,ch in enumerate(bsym)}
	
	rn = reversed(n)
	res = 0
	mult = 1
	for ch in rn:
		res += bsym_inv[ch]*mult
		mult *= len(bsym)
	
	return res

def dec_to_b(n,bsym):
	res = ''
	while n > 0:
		idx = n % len(bsym)
		res = bsym[idx] + res
		n -= idx
		n //= len(bsym)
	return res

'''
p prime modulus
g primitive root mod p
h in Z/pZ
n is number of trials
'''
def bday_attack_disclog(p,g,h,n=None):
	if n is None:
		n = int(np.sqrt(p) * 1.25)
	mr = IntegerModRing(p)
	g = mr(g)
	h = mr(h)
	
	g_exp_try = list(range(p))
	np.random.shuffle(g_exp_try)
	g_exp_try = g_exp_try[:n]
	gh_exp_try = list(range(p))
	np.random.shuffle(gh_exp_try)
	gh_exp_try = gh_exp_try[:n]
	
	
	pg = {(g**exp).x:exp for exp in g_exp_try}
	pgh = {(h*(g**(p-1-exp))).x:exp for exp in gh_exp_try}
	
	g_exp = None
	gh_exp = None
	for g_res,gh_res in zip(pg,pgh):
		if g_res in pgh:
			g_exp = pg[g_res]
			gh_exp = pgh[g_res]
			break
		if gh_res in pg:
			g_exp = pg[gh_res]
			gh_exp = pgh[gh_res]
			break
	
	if g_exp is not None:
		return ModInteger(g_exp,p-1)+gh_exp

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

def nums_to_ascii(nums,block_size):
	msg = ''
	i_block_width = 3  # technically floor(log10(256))+1 but that's just 3
	s_block_width = (int(math.log10(block_size)) + 1) // i_block_width
	for num in nums:
		#pad with zeroes
		numst = str(num)
		while (len(numst) % i_block_width) != 0:
			numst = '0' + numst
		for chsq in [numst[i:i+i_block_width] for i in range(0,len(numst),i_block_width)]:
			numv = int(chsq)
			msg += chr(numv)
	return msg

	
def el_gamal_enc(p,g,msg,pub_key):
	mr = IntegerModRing(p)
	g = mr(g)
	
	#convert msg into numbers
	m_nums = ascii_to_nums(msg,p)
	
	#do all the mults
	rts = []
	for m_num in m_nums:
		#generate a k
		k = np.random.randint(0,p)
		r = g**k
		t = m_num*(mr(pub_key))**(k)
		
		rts.append((r.x,t.x))
	
	return rts


'''
m = t*r^(-a)
'''
def el_gamal_dec(p,rts,a,print_chvals=False):
	mr = IntegerModRing(p)
	rts = [(mr(r),mr(t)) for r,t in rts]
	
	ma = p-1-a
	nums = []
	for r,t in rts:
		m = t*r**(ma)
		#make this into a string, zero-pad, and convert to ascii
		nums.append(m.x)
	return nums_to_ascii(nums,p)
	

def freq_hist(dat,base_alphabet=None,normed=False):
	if base_alphabet is None:
		hist = {}
	else:
		hist = {base:0 for base in base_alphabet}
	for item in dat:
		if item in hist:
			hist[item] += 1
		else:
			hist.update({item:1})
	
	if normed:
		for k in hist:
			hist[k] /= len(dat)
	
	if base_alphabet is None:
		return hist
	
	ahist = []
	for ch in base_alphabet:
		ahist.append(hist[ch])
	
	return ahist


def vig_enc(msg,key):
	msg = remove_non_alphanum(msg)
	key = remove_non_alphanum(key)
	ctxt = ''
	ki = 0
	for ch in msg:
		cch = ((ord(ch)-65) + (ord(key[ki])-65)) % 26
		cch += 65
		cch = chr(cch)
		ctxt += cch
		ki = (ki + 1) % len(key)
	
	return ctxt


def vig_dec(ctxt,key):
	msg = ''
	ki = 0
	for ch in ctxt:
		mch = modi(((ord(ch)-65) - (ord(key[ki])-65)),26)#could be negative, so use a negative-safe mod function
		mch += 65
		mch = chr(mch)
		msg += mch
		ki = (ki + 1) % len(key)
	
	return msg

def coincidence_count(ctxt,offset):
	#ignore the first offset chars in ctxt
	coins = 0
	for i in range(len(ctxt)-offset):
		ch0 = ctxt[i]
		ch1 = ctxt[i+offset]
		if ch0 == ch1:
			coins += 1
	return coins

def get_likely_keylen(ctxt):
	#try a bunch of offsets (between 1 and len(ctxt)) and figure out which one has the highest coincidence counts
	coins = []
	for offset in range(1,int(np.ceil(len(ctxt)/2))):
		coins.append(coincidence_count(ctxt,offset))
	
	#want to "comb" through the coins array until we get a likely keylen (basically, abuse the fact that multiples of the len will also have spikes)
	prev_val = float('inf')
	for offset in range(1,int(np.ceil(len(ctxt)/2))):
		t_a = coins[offset-1:-1:offset]
		val = sum(t_a)/(offset*len(t_a))
		#for a uniform distrib, this ^^^ should be monotonically decreasing, so pick the first length that is a jump up from the previous one
		if val > prev_val:
			return offset
		prev_val = val
	
	#keep this as a fallback
	return np.argmax(coins) + 1#This is the most problematic part (we really want the first big peak, not the biggest peak)

'''
blind letter distribution of some string
'''
def get_letter_distrib(s):
	hist = [0 for _ in range(26)]
	for ch in s:
		hist[ord(ch)-65] += 1
	return np.array(hist)

'''
positive amount shifts left (A becomes B), negative shifts right (A becomes Z)
'''
def shift_freq_distrib(dist,amt):
	ndist = []
	for val in dist[amt:]:
		ndist.append(val)
	for val in dist[:amt]:
		ndist.append(val)
	return np.array(ndist)

'''
get most likely shift of some string s

do this by dotting a shifted frequency distrib of this string with the base english one (max dotp is most likely shift)
'''
def get_shift_of(s):
	freq_distrib = get_letter_distrib(s)
	best_shift = 0
	best_shift_val = -float('inf')
	for shift in range(0,26):
		shifted_distrib = shift_freq_distrib(freq_distrib,shift)
		shift_val = shifted_distrib @ ENGL_FREQ
		if shift_val > best_shift_val:
			best_shift = shift
			best_shift_val = shift_val
	
	return best_shift

def get_vigenere_key(ctxt,verbose=False):
	#first figure out the keylen
	keylen = get_likely_keylen(ctxt)
	if verbose:
		print('MOST LIKELY KEYLEN: {}'.format(keylen))
	
	#then split the ctxt up according to that
	ctxts = [ctxt[i:-1:keylen] for i in range(keylen)]
	
	#then figure out the key for each of those
	reconstructed_key = ''
	for partial_text in ctxts:
		key_this = get_shift_of(partial_text)
		reconstructed_key += chr(key_this + 65)
	
	if verbose:
		#show the most probable reconstructed message
		print('reconstructed message: {}'.format(vig_dec(ctxt,reconstructed_key)))
	
	
	return reconstructed_key

def remove_non_alphanum(s):
	s = s.upper()
	return ''.join([ch for ch in s if ((ord(ch) < 65+26) and (ord(ch) >= 65))])

def caesar_enc(msg,key):
	msg = remove_non_alphanum(msg)
	key = remove_non_alphanum(key)
	mr = IntegerModRing(26)
	keyi = mr(ord(key) - 65)
	
	ctxt = ''
	for ch in msg:
		nch = keyi + (ord(ch) - 65)
		ctxt += chr(nch.x + 65)
	
	return ctxt

def caesar_dec(ctxt,key):
	key = remove_non_alphanum(key)
	mr = IntegerModRing(26)
	keyi = mr(ord(key) - 65)
	
	msg = ''
	for ch in ctxt:
		nch = (-keyi) + (ord(ch) - 65)
		msg += chr(nch.x + 65)
	
	return msg

def caesar_crack(ctxt):
	#basically the same as with vigenere but we only have to do it once
	key = get_shift_of(ctxt)
	return chr(key+65)


def print_all_caesars(ctxt):
	#just show me all the possibilities
	for key in ENGL_FREQ_DICT:
		print(caesar_dec(ctxt,key))

"""
RSA
"""

def RSA_keygen_dec(prime_digits):
	return RSA_keygen(int(prime_digits*(1/np.log10(2))))


'''
get a random prime that is at least this big and at most this big
'''
def get_prime_factor(width):
	#build a random number that is that wide
	rn = '1'
	for _ in range(width-2):
		rn += np.random.choice(['0','1'])

	rn += '1'
	rn = int(rn,base=2)
	#find the next prime (including rn)
	return next_prime(rn,check_n=True)


'''
figure out a large prime factor, then multiply it by another integer and add 1 and see if we get a prime
'''
def make_large_prime(prime_bits,exclude=0,max_iterations=100):
	for _ in range(max_iterations):
		# pick random prime greater than root(p-1)
		phi_pm1_large_fac = get_prime_factor(prime_bits-1)
		for k in PRIMES_1K[:5]:
			p = (phi_pm1_large_fac*k) + 1
			if (p != exclude) and isprime(p):
				return p,phi_pm1_large_fac,k
	raise AttributeError('{} iterations insufficient'.format(max_iterations))


def RSA_keygen(prime_bits):
	#now make p and q
	p,pm1_large_fac,pm1_small_fac = make_large_prime(prime_bits)
	q,qm1_large_fac,qm1_small_fac = make_large_prime(prime_bits,exclude=p)
	#pick random d (invertible mod (p-1)(q-1))
	phi_n = (p-1)*(q-1)
	d = random.randint(2,(p-1)*(q-1))
	while not gcd(phi_n,d) == 1:
		d = random.randint(2, (p - 1) * (q - 1))
	#e is its inverse
	e = ModInteger(1,phi_n,n_is_prime=False)/d
	print('p={}'.format(p))
	print('q={}'.format(q))
	return (d,(p*q,e.x))


'''
both encrypt and decrypt use this

msg is a list of integers
'''
def RSA_en_decrypt(msg,de,n):
	mr = IntegerModRing(n,n_is_prime=False)
	nums = [(mr(x)**de).x for x in msg]#actual de/encryption
	return nums

'''
dealing with string stuff necessitates this
'''
def RSA_encrypt(msg,e,n):
	nums = ascii_to_nums(msg,n)
	nums = RSA_en_decrypt(nums,e,n)
	return nums

def RSA_decrypt(nums,d,n):
	nums = RSA_en_decrypt(nums,d,n)
	return nums_to_ascii(nums,n)

"""
END RSA
"""

'''
p-minus-one factoring method (in class 2020-10-09)
'''
def pm1_factor(n,B_init=None,max_iterations=100,B_inc=None):
	if B_init is None:
		B_init = int(n**(1/6))+1
	if B_inc is None:
		B_inc = int(((n**(1/6))+1)/8)
	mr = IntegerModRing(n)
	#similar to MR
	#pick a in Z/nZ>1
	a = mr(2)#maybe pick something else? wiki says it doesn't matter as long as n is odd (if it isn't it's probably not that hard to factor)
	B = B_init
	for _ in range(max_iterations):
		#compute a chain mod n of
			#a -sq-> a^2 -qbe-> a^(2*3) -...-> a^(B!) =: b
		b = a*a
		for exp in range(2,B):
			b = b**exp
		#try d = gcd(b-1,n) (check if this is nontrivial, repeat if it isn't
		d = gcd(b.x-1,n)
		if d == n:
			if B <= 2:
				return None
			B -= B_inc
		elif d == 1:
			if B >= n:
				return None
			B += B_inc
		else:
			return d,n//d

'''
domain-specific attack
'''
def ten_dig_RSA_attack(n,digits):
	#start with the first possible (k digit) prime
	p = next_prime(10**(digits-1)+1,check_n=True)
	bound = int(10**(digits)) >> 1
	while ((n % p) != 0) and (p < bound):
		#increase p and try again
		p = next_prime(p)
	return p,(n//p)

def trial_division_RSA(n):
	#start with the first possible (k digit) prime
	p = int(math.sqrt(n)+1)
	while (p > 0) and ((n % p) != 0):
		#increase p and try again
		p -= 2
	return p,(n//p)

'''
sounds great, doesn't work
'''
def bday_attack_RSA(c,n,e,list_size=None,B=256,root_coef=1.25):
	rB = int(np.sqrt(B)*root_coef)

	mr = IntegerModRing(n)
	c = mr(c)

	if list_size is None:
		#just do all up to root B
		xs = {(c * x ** (-e)).x:x for x in range(2,rB)}
		ys = {(mr(y) ** (e)).x:y for y in range(2,rB)}
	else:
		xs = {(c * x ** (-e)).x for x in [random.randint(2, rB) for _ in range(list_size)]}
		ys = {(mr(y) ** (e)).x:y for y in [random.randint(2, rB) for _ in range(list_size)]}

	for xr in xs:
		if xr in ys:
			#collision found
			#c cong (xy)^e mod n
			x = xs[xr]
			y = ys[xr]
			return (mr(x)*y)**e