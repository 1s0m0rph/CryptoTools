from ModularArith import *

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
	while not ext_eucl_int(phi_n,d,gcd_only=True) == 1:
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
	return nums_to_ascii(nums)