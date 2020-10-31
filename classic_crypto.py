from crypto_tools.ModularArith import *

def vig_enc(msg, key):
	msg = remove_non_alphanum(msg)
	key = remove_non_alphanum(key)
	ctxt = ''
	ki = 0
	for ch in msg:
		cch = ((ord(ch)-65)+(ord(key[ki])-65))%26
		cch += 65
		cch = chr(cch)
		ctxt += cch
		ki = (ki+1)%len(key)

	return ctxt


def vig_dec(ctxt, key):
	msg = ''
	ki = 0
	for ch in ctxt:
		mch = modi(((ord(ch)-65)-(ord(key[ki])-65)), 26)  #could be negative, so use a negative-safe mod function
		mch += 65
		mch = chr(mch)
		msg += mch
		ki = (ki+1)%len(key)

	return msg


def coincidence_count(ctxt, offset):
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
	for offset in range(1, int(np.ceil(len(ctxt)/2))):
		coins.append(coincidence_count(ctxt, offset))

	#want to "comb" through the coins array until we get a likely keylen (basically, abuse the fact that multiples of the len will also have spikes)
	prev_val = float('inf')
	for offset in range(1, int(np.ceil(len(ctxt)/2))):
		t_a = coins[offset-1:-1:offset]
		val = sum(t_a)/(offset*len(t_a))
		#for a uniform distrib, this ^^^ should be monotonically decreasing, so pick the first length that is a jump up from the previous one
		if val > prev_val:
			return offset
		prev_val = val

	#keep this as a fallback
	return np.argmax(
		coins)+1  #This is the most problematic part (we really want the first big peak, not the biggest peak)


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


def shift_freq_distrib(dist, amt):
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
	for shift in range(0, 26):
		shifted_distrib = shift_freq_distrib(freq_distrib, shift)
		shift_val = shifted_distrib@ENGL_FREQ
		if shift_val > best_shift_val:
			best_shift = shift
			best_shift_val = shift_val

	return best_shift


def get_vigenere_key(ctxt, verbose=False):
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
		reconstructed_key += chr(key_this+65)

	if verbose:
		#show the most probable reconstructed message
		print('reconstructed message: {}'.format(vig_dec(ctxt, reconstructed_key)))

	return reconstructed_key


def caesar_enc(msg, key):
	msg = remove_non_alphanum(msg)
	key = remove_non_alphanum(key)
	mr = IntegerModRing(26)
	keyi = mr(ord(key)-65)

	ctxt = ''
	for ch in msg:
		nch = keyi+(ord(ch)-65)
		ctxt += chr(nch.x+65)

	return ctxt


def caesar_dec(ctxt, key):
	key = remove_non_alphanum(key)
	mr = IntegerModRing(26)
	keyi = mr(ord(key)-65)

	msg = ''
	for ch in ctxt:
		nch = (-keyi)+(ord(ch)-65)
		msg += chr(nch.x+65)

	return msg


def caesar_crack(ctxt):
	#basically the same as with vigenere but we only have to do it once
	key = get_shift_of(ctxt)
	return chr(key+65)


def print_all_caesars(ctxt):
	#just show me all the possibilities
	for key in ENGL_FREQ_DICT:
		print(caesar_dec(ctxt, key))