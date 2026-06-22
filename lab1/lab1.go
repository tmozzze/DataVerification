package lab1

import (
	"math"
	"math/cmplx"
)

type DataVerifier struct{}

func gcd(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}

	return a
}

func (dv DataVerifier) EulerByDef(n int) int {
	if n < 1 {
		return 0
	}

	res := 0
	for i := 1; i <= n; i++ {
		if gcd(i, n) == 1 {
			res++
		}
	}
	return res
}

func (dv DataVerifier) EulerByFactorization(n int) int {
	if n < 1 {
		return 0
	}
	if n == 1 {
		return 1
	}

	// factorization
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

func (dv DataVerifier) EulerByFourier(n int) int {
	sum := .0
	for k := 1; k <= n; k++ {
		angle := 2.0 * math.Pi * float64(k) / float64(n)

		term := float64(gcd(k, n)) * math.Cos(angle)
		sum += term
	}

	return int(math.Round(sum))
}

func (dv DataVerifier) RootsOfUnity(n int) []complex128 {
	roots := make([]complex128, n)

	for j := 0; j < n; j++ {
		angle := 2.0 * math.Pi * float64(j) / float64(n)

		roots[j] = complex(math.Cos(angle), math.Sin(angle))
	}
	return roots
}

func (dv DataVerifier) PrimitiveRootsOfUnity(n int) []complex128 {
	var primitives []complex128
	for j := 0; j < n; j++ {
		if gcd(j, n) == 1 {
			angle := 2.0 * math.Pi * float64(j) / float64(n)
			primitives = append(primitives, complex(math.Cos(angle), math.Sin(angle)))
		}
	}
	return primitives
}

func (dv DataVerifier) VandermondeMatrices(n int, primitiveRoot complex128) ([][]complex128, [][]complex128) {
	V := make([][]complex128, n)
	Vinv := make([][]complex128, n)

	for j := 0; j < n; j++ {
		V[j] = make([]complex128, n)
		Vinv[j] = make([]complex128, n)
		for k := 0; k < n; k++ {

			power := float64(j * k)

			V[j][k] = cmplx.Pow(primitiveRoot, complex(power, 0))

			Vinv[j][k] = cmplx.Pow(primitiveRoot, complex(-power, 0)) / complex(float64(n), 0)
		}
	}
	return V, Vinv
}

func (dv DataVerifier) FurierTransform(matrix [][]complex128, vector []complex128) []complex128 {
	n := len(vector)
	result := make([]complex128, n)

	for i := 0; i < n; i++ {
		var sum complex128 = 0
		for j := 0; j < n; j++ { // Идем по элементам строки и вектора
			sum += matrix[i][j] * vector[j] // Умножаем и складываем
		}
		result[i] = sum // Записываем результат в новый массив
	}
	return result
}
