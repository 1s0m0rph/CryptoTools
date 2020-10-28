from crypto_tools.ModularArith import *

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
		d = ext_eucl_int(b.x-1,n,gcd_only=True)
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

'''
not just the factor base, but we need it to have square roots for n modulo p
also tell us what those roots are

source for quadradic residues: Weissman pp 199
'''
def get_qs_factor_base(n,B,kmax):
	ptb = primes(B)[1:]#2 sucks so we're not going to test it (it works as long as n is odd anyway)
	if (n & 1) == 1:#2 works
		fb = {(2,1):{ModInteger(1,2)}}
		order = [(2,1)]  #what order do we go in when iterating through this
		#2 will have to be done manually (i.e. no powers)
	else:#2 doesn't work
		fb = {}
		order = []
	for p in ptb:
		#make sure n has a square root modulo p
		#use euler criterion:
		#a is square mod p iff a^((p-1)/2) == 1 mod p
		nmp = ModInteger(n,p)
		if (p % 4) == 3:#technically a partial duplication of ModInteger.sqrt but it is being used differently here
			sqrt = nmp**((p+1)//4)
			if sqrt**2 == nmp:
				fb.update({(p,1):{sqrt,-sqrt}})
		elif (p % 8) == 5:
			#here there are actually (as many as) 4 sqrts
			sqrt1 = nmp**((n+3)//8)
			sqrt2 = nmp**((n+3)//8)*ModInteger(2,p)**((n-1)//4)
			if sqrt1**2 == nmp:
				fb.update({(p,1):{sqrt1,-sqrt1,sqrt2,-sqrt2}})
		#else it's 1 mod 8 and no general form solution exists, so let's just pretend this number doesn't exist

		if (p,1) in fb:
			#now modify the sieve so that we don't have to do any trial division
			#this means we need to add powers of the primes in, which can be done using hensel's lemma
			#we'll have solutions of the form s = r - (r^2 - n)*a mod p^k+1 where a == (2r)^-1 mod p
			#we only have to do this until p^k > kmax
			k = 2
			npki = p**2
			inv_mod_p = {}#just for efficiency's sake (caching)
			order.append((p,1))
			while npki <= kmax:
				order.append((p,k))
				fb.update({(p,k):set()})
				for sqrt in fb[(p,k-1)]:
					#add on a root mod p^k = npki
					#first find a
					ntinv = ModInteger(sqrt.x,p)*2
					if ntinv.x in inv_mod_p:
						a = inv_mod_p[ntinv.x]
					else:
						a = (ntinv**(-1)).x
						inv_mod_p.update({ntinv.x:a})

					#now find s
					s = ModInteger(sqrt,npki) - ModInteger(a,npki)*(sqrt.x**2 - n)
					fb[(p,k)].add(s)
				npki *= p
				k += 1

	return fb,order


'''
quadratic sieve implementation
'''
def quad_sieve(n,B,kmin=None,knum=None):
	if kmin is None:
		kmin = int(np.ceil(np.sqrt(n)))
	if knum is None:
		knum = B
	kmax = kmin + knum
	fb,order = get_qs_factor_base(n,B,kmax)#includes primes with roots as well as their powers up to kmin+knum
	facts = {k:None for k in range(kmin,kmax+1)}#maps k onto factorization if and only if k^2 is B-smooth with this factor base
	for i,(p,e) in enumerate(order):
		pe = p**e
		k_start = kmin//pe
		offset = pe*k_start#0 mod p^e
		while offset < kmax:
			for sqrt in fb[(p,e)]:
				nkeep = offset + sqrt.x
				target = nkeep**2 - n
				if (nkeep in facts) and (target >= pe):
					if facts[nkeep] is None:
						facts[nkeep] = [{p:0 for _ in range(p)},1,False]#(exponents start all at 0,product so far,product at target)
					if not facts[nkeep][2]:
						if p not in facts[nkeep][0]:
							facts[nkeep][0].update({p:0})
						if p == 2:#2 has to go backwards because it is evil
							k,rem = get_int_2km_fac(target)
							facts[nkeep][0][2] = k
							facts[nkeep][1] <<= k
							if rem == 1:
								facts[nkeep][2] = True
						elif nkeep >= kmin:
							facts[nkeep][0][p] += 1
							facts[nkeep][1] *= p
							if facts[nkeep][1] == target:
								facts[nkeep][2] = True
			offset += pe

	facts = [[k,facts[k][0]] for k in facts if (facts[k] is not None) and (facts[k][2])]

	#now that we have these, we can ACTUALLY factor n (probably)
	#we need to find some combination of these vectors (i.e. the exponent vectors) that sums to 0 mod 2
	seen_ones = {}
	ps = list({p for p,_ in fb})
	A = Matrix([[ModInteger(facts[i][1][p],2) if p in facts[i][1] else ModInteger(0,2) for i in range(len(facts))] for p in ps],dtype=int)
	Z = A.nullspace()#finds many solutions simultaneously -- just pick the first one
	if len(Z) == 0:
		raise AttributeError("No solutions found, try increasing B")
	sol = set()
	for i,x in enumerate(Z[0]):
		if ModInteger(int(x),2) == 1:
			sol.add(i)



	#now we should have some y such that y^2 == x^2 mod n
	#find x and y
	x = ModInteger(1,n)
	yexps = {p:0 for p in ps}
	for i in sol:
		x *= facts[i][0]#x multiplies by k
		for p in facts[i][1]:
			yexps[p] += facts[i][1][p]
	y = ModInteger(1,n)
	for p in ps:
		y *= ModInteger(p,n)**(yexps[p]//2)

	fac = ext_eucl_int((x-y).x,n,gcd_only=True)
	return fac,n//fac