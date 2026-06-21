package lab1

type DataVerifier struct{}

func gcd(a, b int) int {
	for b != 0 {
		a, b = b, a%b
	}

	return a
}

func (dv DataVerifier) EulerByDef(n int) int {
	res := 0
	for i := 1; i <= n; i++ {
		if gcd(i, n) == 1 {
			res++
		}
	}
	return res
}
