"""
Lattice-based cryptography. Mostly just ring- learning with errors
"""

from crypto_tools.FiniteFields import *


class RLWE:
	"""
	Used to set up an RLWE-system with the prior public info
	"""

	def __init__(self,p:int,n:int,k:int,smallness_bound:int,**ffp_kwargs):
		self.p = p
		self.n = n
		self.msg_lim = p//k#p is prime so highest possible number is just floor(p/k) (unless k is 1, which wouldn't work anyway)
		self.k = k
		self.m = FiniteFieldPoly(p,[1]+([0]*(n-1))+[1],**ffp_kwargs)#assume the default modulus of x^n + 1
		self.R = smallness_bound
		self.ffp_kwargs = ffp_kwargs

		self.FFPM = FiniteFieldModM(p,self.m,**ffp_kwargs)

	'''
	allow for predefined private keys
	'''
	def generate_keypair(self,private_key:Union[FiniteFieldPoly,FiniteFieldPolyModM]=None) -> (FiniteFieldPolyModM,FiniteFieldPolyModM,FiniteFieldPolyModM):
		if private_key is None:
			private_key = FFP_with_random_coefs(self.p,self.n,mag_bound=self.R,m=self.m,**self.ffp_kwargs)
		if type(private_key) != FiniteFieldPolyModM:
			private_key = self.FFPM(private_key)
		elif private_key.m != self.m:
			raise AttributeError("Polynomial modulus of private key is not the same as system modulus.")
		assert(private_key.p == self.p)
		assert(private_key.poly.dgr <= self.n)

		#now set up the public key (a,b) where b = as + e
		a = FFP_with_random_coefs(self.p,self.n-1,m=self.m,**self.ffp_kwargs)#random a (NOT small)
		e = FFP_with_random_coefs(self.p,self.n-1,mag_bound=self.R,m=self.m,**self.ffp_kwargs)#random error (small)

		#now create b
		b = (a*private_key) + e

		return private_key,a,b


	def check_polys_are_in_system(self,polys:List[FiniteFieldPolyModM]):
		for x in polys:
			try:
				assert (x.p == self.p)
				assert (x.m == self.m)
				assert (x.poly.dgr <= self.n)
			except AssertionError:
				raise AttributeError("Polynomial {} is not in the RLWE system (one of its attributes faild to match)".format(x))

	def encrypt_single_message_to_key(self,a:FiniteFieldPolyModM,b:FiniteFieldPolyModM,m:int):
		if (m < 0) or (m > self.msg_lim):
			raise AttributeError("Message {} out of bounds for p={}, k={} (max message size {})".format(m,self.p,self.k,self.msg_lim))


		self.check_polys_are_in_system([a,b])

		#pick ephemeral key r (small)
		r = FFP_with_random_coefs(self.p,self.n-1,self.R,self.m,**self.ffp_kwargs)
		#pick two small errors
		e1 = FFP_with_random_coefs(self.p,self.n-1,self.R,self.m,**self.ffp_kwargs)
		e2 = FFP_with_random_coefs(self.p,self.n-1,self.R,self.m,**self.ffp_kwargs)

		#define v and w (ciphertext)
		v = (a*r) + e1
		w = (b*r) + e2 + (m*self.k)

		return v,w

	def encrypt(self,a:FiniteFieldPolyModM,b:FiniteFieldPolyModM,msg:str):
		max_strlen = int(math.floor(math.log10(self.msg_lim-1))+1)
		m_nums = ascii_to_nums_nomod(msg,max_strlen//3)
		ctxts = []

		for m in m_nums:
			ctxt = self.encrypt_single_message_to_key(a,b,m)
			ctxts.append(ctxt)

		return ctxts


	def decrypt_single_message(self,v:FiniteFieldPolyModM,w:FiniteFieldPolyModM,private_key:FiniteFieldPolyModM):
		#check all the things
		self.check_polys_are_in_system([w,v,private_key])

		#start by taking the polynomial (w - vs)
		unrounded = w - (v*private_key)
		unround_const = unrounded.poly[0].x#sets all the higher-degree (i.e. nonconstant terms) to zero
		m = int(np.round(unround_const / self.k))#just use numpy to round instead of doing it ourselves and messing it up

		return m

	def decrypt(self,wvs:List[Union[Tuple[FiniteFieldPolyModM],List[FiniteFieldPolyModM]]],private_key:FiniteFieldPolyModM) -> str:
		ms = []

		for w,v in wvs:
			m = self.decrypt_single_message(w,v,private_key)
			ms.append(m)

		msg = nums_to_ascii(ms)
		return msg