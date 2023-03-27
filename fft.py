import math
import utils
from utils import random_input
from utils import complex_sum
from utils import complex_minus
from utils import complex_product
from utils import complex_l1_dis
from utils import arr_sum
from utils import arr_minus
from utils import arr_scale

# complex numbers are represented as tuple a + bi -> (a, b)


w_twiddle = dict({})
w_twiddle_conj = dict({})
ws_twiddle = dict({})
ws_twiddle_conj = dict({})
t_twiddle = dict({})
t_twiddle_conj = dict({})

# although these are always real, we store them as a tuple to be consistent and just say they have imaginary part 0
s2_const = dict({})
s4_const = dict({})

def init_constants():
	def const_s(N, k):
		# assuming N is a power of 2
		if N < 4:
			return 1

		k4 = k % (N / 4)
		angle = 2 * math.pi * k4 / N
		rec = const_s(N/4, k4)
		if k4 <= N/8:
			return math.cos(angle) * rec
		return math.sin(angle) * rec

	for n in range(int(math.log2(N)) + 1):
		temp_N = 2 ** n
		for k in range(temp_N + 1):
			real_part = math.cos(2 * math.pi * k / temp_N)
			imag_part = math.sin(-2 * math.pi * k / temp_N)

			w_twiddle[(temp_N, k)] =      (real_part,      imag_part)
			w_twiddle_conj[(temp_N, k)] = (real_part, -1 * imag_part)

			s_N_k = const_s(temp_N, k)
			s_N4_k = const_s(temp_N / 4, k)

			ws_twiddle[(temp_N, k)] =      (real_part * s_N4_k,      imag_part * s_N4_k)
			ws_twiddle_conj[(temp_N, k)] = (real_part * s_N4_k, -1 * imag_part * s_N4_k)

			t_twiddle[(temp_N, k)] =      (real_part * s_N4_k / s_N_k,      imag_part * s_N4_k / s_N_k)
			t_twiddle_conj[(temp_N, k)] = (real_part * s_N4_k / s_N_k, -1 * imag_part * s_N4_k / s_N_k)

			if temp_N < N/2:
				s2_const[(temp_N, k)] = (const_s(temp_N, k) / const_s(2 * temp_N, k), 0)
			if temp_N < N/4:
				s4_const[(temp_N, k)] = (const_s(temp_N, k) / const_s(4 * temp_N, k), 0)

class split_radix:

	def compute_FFT(nums):
		return split_radix.compute(nums)

	def compute(nums):
		local_N = len(nums)
		if local_N == 2:
			# 2x2 fft base case, which is also 2x2 hadamard ase case
			return [complex_sum(nums[0], nums[1]), complex_minus(nums[0], nums[1])]
		elif local_N == 1:
			return nums


		u = split_radix.compute(nums[0 : local_N // 2])
		zk = split_radix.compute(nums[local_N // 2 : 3 * local_N // 4])
		zpk = split_radix.compute(nums[3 * local_N // 4 : local_N])

		'''
		print(u)
		print(zk)
		print(zpk)
		'''

		y = [(0, 0)] * local_N

		for k_0 in range(local_N // 4):
			k_1 = k_0 + local_N // 4
			k_2 = k_1 + local_N // 4
			k_3 = k_2 + local_N // 4

			indices = [k_0, k_1, k_2, k_3]

			# multiply zk and zpk in twiddle first, in place

			zk[k_0] = complex_product(zk[k_0], w_twiddle[(local_N, k_0)])
			zpk[k_0] = complex_product(zpk[k_0], w_twiddle_conj[(local_N, k_0)])

			z_sum = complex_sum(zk[k_0], zpk[k_0])
			z_dif = complex_minus(zk[k_0], zpk[k_0])

			for i in range(len(indices)):
				y[indices[i]] = u[indices[i % 2]]

			'''
			print(y)
			print(z_sum)
			print(z_dif)
			'''

			y[k_0] = complex_sum(y[k_0], z_sum)
			y[k_2] = complex_minus(y[k_2], z_sum)

			z_diff_times_negi = complex_product(z_dif, (0, -1)) # beware when counting ops

			y[k_1] = complex_sum(y[k_1], z_diff_times_negi)
			y[k_3] = complex_minus(y[k_3], z_diff_times_negi)

		return y

class johnson_frigo:

	def compute_FFT(nums):
		return johnson_frigo.compute(0, nums)

	def compute(func_num, nums):
		local_N = len(nums)
		# computes the johnson frigo FFT of nums assuming nums has power of 2 length
		# func_num = 0 is F, = 1 is FS, = 2 is FS2, = 3 is FS4
		if local_N == 2:
			# 2x2 fft base case, which is also 2x2 hadamard ase case
			(nums_0_real, nums_0_imag) = nums[0]
			(nums_1_real, nums_1_imag) = nums[1]
			ans = [complex_sum(nums[0], nums[1]), complex_minus(nums[0], nums[1])]

			if func_num == 3:
				ans[1] = complex_product(ans[1], s4_const[(local_N, 1)])

			return ans

		elif local_N == 1:
			return nums

		routing = [0, 2, 3, 2] # routing[func_num] decides which func_num to use next

		u = johnson_frigo.compute(routing[func_num], nums[0 : local_N // 2])
		zk = johnson_frigo.compute(1, nums[local_N // 2 : 3 * local_N // 4])
		zpk = johnson_frigo.compute(1, nums[3 * local_N // 4 : local_N])

		'''
		print(u)
		print(zk)
		print(zpk)
		'''

		y = [(0, 0)] * local_N

		for k_0 in range(local_N // 4):
			k_1 = k_0 + local_N // 4
			k_2 = k_1 + local_N // 4
			k_3 = k_2 + local_N // 4

			indices = [k_0, k_1, k_2, k_3]

			if func_num == 0:
				zk[k_0] = complex_product(zk[k_0], ws_twiddle[(local_N, k_0)])
				zpk[k_0] = complex_product(zpk[k_0], ws_twiddle_conj[(local_N, k_0)])

			else:
				zk[k_0] = complex_product(zk[k_0], t_twiddle[(local_N, k_0)])
				zpk[k_0] = complex_product(zpk[k_0], t_twiddle_conj[(local_N, k_0)])

			'''
			print("after scale")
			print(zk)
			print(zpk)
			'''

			z_sum = complex_sum(zk[k_0], zpk[k_0])
			z_dif = complex_minus(zk[k_0], zpk[k_0])

			if func_num == 2:
				z_sum = complex_product(z_sum, s2_const[(local_N, k_0)])
				z_dif = complex_product(z_dif, s2_const[(local_N, k_1)])

			for i in range(len(indices)):
				y[indices[i]] = u[indices[i % 2]]

			'''
			print(y)
			print(z_sum)
			print(z_dif)
			'''

			y[k_0] = complex_sum(y[k_0], z_sum)
			y[k_2] = complex_minus(y[k_2], z_sum)

			z_diff_times_negi = complex_product(z_dif, (0, -1)) # beware when counting ops

			y[k_1] = complex_sum(y[k_1], z_diff_times_negi)
			y[k_3] = complex_minus(y[k_3], z_diff_times_negi)

			if func_num == 3:
				for i in range(4):
					y[indices[i]] = complex_product(y[indices[i]], s4_const[(local_N, indices[i])])

		return y
		# TODO validation test jf fft

class hadamard:

	def naive_hadamard(nums):
		if len(nums) == 1:
			return nums
		elif len(nums) == 2:
			temp = nums[0]
			nums[0] = complex_sum(temp, nums[1])
			nums[1] = complex_minus(temp, nums[1])
			return nums

		half = len(nums) // 2
		left = hadamard.naive_hadamard(nums[0 : half])
		right = hadamard.naive_hadamard(nums[half : len(nums)])

		for i in range(half):
			nums[i] = complex_sum(left[i], right[i])
			nums[i + half] = complex_minus(left[i], right[i])

		return nums

	def fast_WHT(nums, scale):
		if len(nums) <= 4:
			return arr_scale(hadamard.naive_hadamard(nums), 2 ** scale)

		step = len(nums) // 8

		a = hadamard.fast_WHT(nums[0 : step], scale)
		b = hadamard.fast_WHT(nums[step * 1 : step * 2], scale + 1)
		c = hadamard.fast_WHT(nums[step * 2 : step * 3], scale + 1)
		d = hadamard.fast_WHT(nums[step * 3 : step * 4], scale + 1)
		e = hadamard.fast_WHT(nums[step * 4 : step * 5], scale + 1)
		f = hadamard.fast_WHT(nums[step * 5 : step * 6], scale + 1)
		g = hadamard.fast_WHT(nums[step * 6 : step * 7], scale + 1)
		h = hadamard.fast_WHT(nums[step * 7 : step * 8], scale + 1)

		B_1 = arr_sum(b, c)
		B_2 = arr_sum(d, h)
		B_3 = arr_sum(f, g)

		tot = arr_sum(arr_sum(arr_sum(B_1, B_2), B_3), e)
		tot = arr_scale(tot, 0.5)

		diff = arr_minus(a, tot)

		D = arr_sum(diff, d)
		E = arr_sum(diff, e)
		H = arr_sum(diff, h)

		y_1 = arr_sum(a, tot)
		y_2 = arr_sum(arr_sum(E, c), g)
		y_3 = arr_sum(arr_sum(E, b), f)
		y_4 = arr_sum(E, B_2)
		y_5 = arr_sum(D, B_1)
		y_6 = arr_sum(arr_sum(H, c), f)
		y_7 = arr_sum(arr_sum(H, b), g)
		y_8 = arr_sum(D, B_3)

		return y_1 + y_2 + y_3 + y_4 + y_5 + y_6 + y_7 + y_8

class WHUFFT:

	def H_prime(nums):
		def arr_at_indices(arr, indices):
			# assumes indices are in sorted order
			ans = []
			for index in indices:
				ans.append(arr[index])
			return ans
		partition = WHUFFT.Partition(len(nums))

		y = [(0, 0)] * len(nums)

		for subset in partition:
			selected_input = arr_at_indices(nums, subset)
			selected_output = hadamard.fast_WHT(selected_input, 0)
			for i in range(len(subset)):
				y[subset[i]] = selected_output[i]

		return y

	def Partition(N):
		if N == 1:
			return [[0]]
		elif N == 2:
			return [[0], [1]]

		ans = []
		P1 = WHUFFT.Partition(N // 2)
		P2 = WHUFFT.Partition(N // 4)

		for subset in P1:
			ans.append(subset)

		for subset in P2:
			temp = []
			for index in subset:
				temp.append(index + N // 2)
				temp.append(index + 3 * N // 4)

			temp.sort()

			ans.append(temp)

		return ans

	def compute_FFT(nums):
		return WHUFFT.compute(0, WHUFFT.H_prime(nums))

	def compute(func_num, nums):
		local_N = len(nums)

		if local_N == 2:
			# 2x2 fft base case, which is also 2x2 hadamard ase case
			(nums_0_real, nums_0_imag) = nums[0]
			(nums_1_real, nums_1_imag) = nums[1]
			ans = [complex_sum(nums[0], nums[1]), complex_minus(nums[0], nums[1])]

			if func_num == 3:
				ans[1] = complex_product(ans[1], s4_const[(local_N, 1)])

			return ans

		elif local_N == 1:
			return nums

		routing = [0, 2, 3, 2] # routing[func_num] decides which func_num to use next

		u = WHUFFT.compute(routing[func_num], nums[0 : local_N // 2])
		zk = WHUFFT.compute(1, nums[local_N // 2 : 3 * local_N // 4])
		zpk = WHUFFT.compute(1, nums[3 * local_N // 4 : local_N])

		y = [(0, 0)] * local_N

		for k_0 in range(local_N // 4):
			k_1 = k_0 + local_N // 4
			k_2 = k_1 + local_N // 4
			k_3 = k_2 + local_N // 4

			indices = [k_0, k_1, k_2, k_3]

			# declare variables r, r'
			real_twiddle = 0
			imag_twiddle = 0

			if func_num == 0:
				(real_twiddle, imag_twiddle) = ws_twiddle[(local_N, k_0)]
			else:
				(real_twiddle, imag_twiddle) = t_twiddle[(local_N, k_0)]

			(b_plus_c, bprime_plus_cprime) = zk[k_0]
			(b_minus_c, bprime_minus_cprime) = zpk[k_0]

			D = real_twiddle * b_plus_c + imag_twiddle * -1 * bprime_minus_cprime;
			E = real_twiddle * bprime_plus_cprime + imag_twiddle * b_minus_c;
			F = imag_twiddle * b_plus_c + real_twiddle * bprime_minus_cprime;
			G = imag_twiddle * bprime_plus_cprime + real_twiddle * -1 * b_minus_c;

			if func_num == 2:
				D *= s2_const[(local_N, k_0)][0] # real part of the tuple
				E *= s2_const[(local_N, k_0)][0] # real part of the tuple
				F *= s2_const[(local_N, k_1)][0] # real part of the tuple
				G *= s2_const[(local_N, k_1)][0] # real part of the tuple

			a_real = u[k_0][0]
			a_imag = u[k_0][1]
			z_real = u[k_1][0]
			z_imag = u[k_1][1]

			y[k_0] = (a_real + D, a_imag + E)
			y[k_1] = (z_real + F, z_imag + G)
			y[k_2] = (a_real - D, a_imag - E)
			y[k_3] = (z_real - F, z_imag - G)

			if func_num == 3:
				for index in indices:
					y[index] = complex_product(y[index], s4_const[(local_N, index)])

		return y

if __name__ == "__main__":
	global N
	N = 524288 # the input size for everything.

	init_constants()

	sr_input = random_input(N)
	jf_input = sr_input.copy()
	nf_input = sr_input.copy()

	print(len(sr_input))

	sr_output = split_radix.compute_FFT(sr_input)
	jf_output = johnson_frigo.compute_FFT(jf_input)
	nf_output = WHUFFT.compute_FFT(nf_input)

	(tot_dis, max_dis) = complex_l1_dis(sr_output, jf_output)

	print("comparing johnson-frigo to split-radix")
	print("distance between outputs: " + str(tot_dis))
	print("max distance on an index: " + str(max_dis))

	(tot_dis, max_dis) = complex_l1_dis(sr_output, nf_output)

	print("comparing WHUFFT to split-radix")
	print("distance between outputs: " + str(tot_dis))
	print("max distance on an index: " + str(max_dis))




