from EllipticCurves import *

"""
main goal here is to implement at least a rudimentary version of SIDH or (ideally) CSIDH
"""


#TODO finish isogeny crypto implementation
'''
TODOS:

how do we generate an initial EC?
	this needs to have the good properties:
		has p+1 points (where p is the modulus) up to change of coordinates
			use j-invariant here (what is?)
	something from lecture: "for good primes, y^2 = x^3 + 1 will work for this ^^"
		from https://math.stackexchange.com/questions/234837/prove-that-there-are-p1-points-on-the-elliptic-curve-y2-x3-1-over-m:
		"there are p+1 points on the EC y^2 = x^3 + 1 over FF_p where p>3 is a prime st p == 2 mod 3"
		so just pick p == 2 mod 3 and this initial will work
		from the csidh paper (pg. 7), we actually want to pick p == 11 mod 12 [only works for FF_p]
		https://csidh.isogeny.org/csidh-20181118.pdf
		alternatively, as parametrized on pg 13 of the above, we could (read: should) use E0: y^2 = x^3 + x over FF_p
		
how do we create an isogeny from one EC to another?
	this needs to have the good properties:
		one for degree l1 AND l2 (kernel is size l1,l2 -- l1 is small hops, l2 is bigs. 2 per EC [symmetry will help here])
	
	see pg 27 of the original csidh paper
'''