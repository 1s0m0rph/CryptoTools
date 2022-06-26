"""
Basic utils that are commonly needed but not really a part of any one module

these will generally be imported by all modules

basic principle with this file: if it CAN go here (i.e. no missing references) it SHOULD go here
"""


#all project imports go here
from typing import Union, List, Tuple

import numpy as np
from sympy import Matrix
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


def long_divide(a,b,non_neg_rem=False):
	"""
	give integers q,r such that b*q + r = a
	
	if non_neg_rem is set, r will always be nonnegative. otherwise, will return the q,r pair with the smallest (absolute value) remainder
	"""

	q = a//b
	adiffr0 = a - q*b
	adiff0 = abs(adiffr0)
	adiffr1 = adiffr0 - b
	adiff1 = abs(adiffr1)
	if (adiff0 < adiff1) or (non_neg_rem):
		return q,adiffr0
	else:
		return q+1,adiffr1


'''
factorize n as 2^k * m with k maximal
'''
def get_int_2km_fac(n,rw=False):
	nbin = list(reversed(bin(n)[2:]))
	nbwidth = len(nbin)
	k = 0
	while (len(nbin) > 0) and (nbin[0] == '0'):
		k += 1
		del nbin[0]
	#remaining part is m
	nbin = ''.join(reversed(nbin))
	m = int(nbin, base=2)
	if rw:
		return k,m,nbwidth
	else:
		return k,m


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


def b_to_dec(n, bsym):
	#invert the bsymbols
	bsym_inv = {ch: i for i, ch in enumerate(bsym)}

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


def ascii_to_nums_nomod(s,s_block_width):
	i_block_width = 3  #technically floor(log10(256))+1 but that's just 3
	m_nums = []
	for chsq in [s[i:i+s_block_width] for i in range(0, len(s), s_block_width)]:
		nums_this = ''
		for ch in chsq:
			tmp_chord = str(ord(ch))
			tmp_chord = ('0'*(i_block_width-len(tmp_chord)))+tmp_chord
			nums_this += tmp_chord
		m_nums.append(int(nums_this))
	return m_nums

def nums_to_ascii(nums):
	msg = ''
	i_block_width = 3  # technically floor(log10(256))+1 but that's just 3
	for num in nums:
		#pad with zeroes
		numst = str(num)
		while (len(numst) % i_block_width) != 0:
			numst = '0' + numst
		for chsq in [numst[i:i+i_block_width] for i in range(0,len(numst),i_block_width)]:
			numv = int(chsq)
			msg += chr(numv)
	return msg


def freq_hist(dat, base_alphabet=None, normed=False):
	if base_alphabet is None:
		hist = {}
	else:
		hist = {base: 0 for base in base_alphabet}
	for item in dat:
		if item in hist:
			hist[item] += 1
		else:
			hist.update({item: 1})

	if normed:
		for k in hist:
			hist[k] /= len(dat)

	if base_alphabet is None:
		return hist

	ahist = []
	for ch in base_alphabet:
		ahist.append(hist[ch])

	return ahist


def remove_non_alphanum(s):
	s = s.upper()
	return ''.join([ch for ch in s if ((ord(ch) < 65+26) and (ord(ch) >= 65))])

def modi(a,n):
	if a < 0:
		return n*(-(a // n)) + a
	return a % n

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
