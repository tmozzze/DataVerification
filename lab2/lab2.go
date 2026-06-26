package lab2

type DataVerifier struct{}

func gcd(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}
	return a
}

func (dv DataVerifier) euler(n int) int {
	result := n
	for p := 2; p*p <= n; p++ {
		if n%p == 0 {
			for n%p == 0 {
				n /= p
			}
			result -= result / p
		}
	}
	if n > 1 {
		result -= result / n
	}
	return result
}

// степень по модулю
func power(base, exp, mod int) int {
	res := 1
	base %= mod
	for exp > 0 {
		if exp%2 == 1 {
			res = (res * base) % mod
		}
		base = (base * base) % mod
		exp /= 2
	}
	return res
}

// простые множители
func primeFactors(n int) map[int]bool {
	factors := make(map[int]bool)
	for i := 2; i*i <= n; i++ {
		if n%i == 0 {
			factors[i] = true
			for n%i == 0 {
				n /= i
			}
		}
	}
	if n > 1 {
		factors[n] = true
	}
	return factors
}

// маалая теорема Ферма
func modInverse(a, m int) int {
	return power(a, m-2, m)
}

func (dv DataVerifier) hasPrimitiveRoot(n int) bool {
	if n == 2 || n == 4 {
		return true
	}
	if n%2 == 0 {
		n /= 2
	}
	p := 2
	for n%p == 0 {
		n /= p
	} // пропускаем двойку
	p = 3
	for p*p <= n {
		if n%p == 0 {
			for n%p == 0 {
				n /= p
			}
			if n > 1 {
				return false
			} // если есть другой простой множитель
		}
		p += 2
	}
	return true
}

func (dv DataVerifier) FindAllPrimitiveRoots(n int) []int {
	if !dv.hasPrimitiveRoot(n) {
		return nil // корней нет
	}

	phi := dv.euler(n)
	phiFactors := primeFactors(phi)
	var firstRoot int

	// 1 корень
	for g := 2; g < n; g++ {
		if gcd(g, n) != 1 {
			continue
		}
		isPrimitive := true
		for p := range phiFactors {
			if power(g, phi/p, n) == 1 {
				isPrimitive = false
				break
			}
		}
		if isPrimitive {
			firstRoot = g
			break
		}
	}

	if firstRoot == 0 {
		return nil
	}

	// остальные корни из первого
	var allRoots []int
	for i := 1; i <= phi; i++ {
		if gcd(i, phi) == 1 {
			root := power(firstRoot, i, n)
			allRoots = append(allRoots, root)
		}
	}
	return allRoots
}

func (dv DataVerifier) VandermondeMatrixMod(size, n int, omega int) [][]int {
	V := make([][]int, size)
	for i := 0; i < size; i++ {
		V[i] = make([]int, size)
		for j := 0; j < size; j++ {
			V[i][j] = power(omega, i*j, n)
		}
	}
	return V
}

func (dv DataVerifier) SlowFurierTrans(matrix [][]int, vector []int, n int) []int {
	size := len(vector)
	result := make([]int, size)
	for i := 0; i < size; i++ {
		sum := 0
		for j := 0; j < size; j++ {
			term := (matrix[i][j] * vector[j]) % n
			sum = (sum + term) % n
		}
		result[i] = sum
	}
	return result
}

func (dv DataVerifier) FastFurierTrans(vector []int, n int, omega int) []int {
	size := len(vector)
	if size == 1 {
		return vector
	}

	// четные и нечетные
	evens := make([]int, size/2)
	odds := make([]int, size/2)
	for i := 0; i < size/2; i++ {
		evens[i] = vector[2*i]
		odds[i] = vector[2*i+1]
	}

	// для половинок
	fftEvens := dv.FastFurierTrans(evens, n, power(omega, 2, n))
	fftOdds := dv.FastFurierTrans(odds, n, power(omega, 2, n))

	result := make([]int, size)
	currentOmega := 1
	for k := 0; k < size/2; k++ {
		t := (currentOmega * fftOdds[k]) % n
		result[k] = (fftEvens[k] + t) % n
		result[k+size/2] = (fftEvens[k] - t + n) % n // +n чтобы не было отрицательных
		currentOmega = (currentOmega * omega) % n
	}
	return result
}

// обратное
func (dv DataVerifier) InvFastFurierTrans(vector []int, n int, omega int) []int {
	size := len(vector)

	// с обратным корнем
	omegaInv := modInverse(omega, n)
	result := dv.FastFurierTrans(vector, n, omegaInv)

	// результат делим на n
	sizeInv := modInverse(size, n)
	for i := 0; i < size; i++ {
		result[i] = (result[i] * sizeInv) % n
	}
	return result
}
