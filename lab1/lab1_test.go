package lab1

import (
	"math/cmplx"
	"testing"
)

const epsilon = 1e-9

func isCopmlexEqual(a, b complex128) bool {
	return cmplx.Abs(a-b) < epsilon
}

func TestEulerFunctions(t *testing.T) {
	dv := DataVerifier{}

	tests := []struct {
		name     string
		n        int
		expected int
	}{
		{"Один", 1, 1},
		{"Простое", 5, 4},
		{"Составное 6", 6, 2},
		{"Степень простго 9", 9, 6},
		{"Составное 10", 10, 4},
		{"Составное число 12", 12, 4},
		{"Составное число 15", 15, 8},
		{"Большое составное 36", 36, 12},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {

			resDef := dv.EulerByDef(tc.n)
			resFac := dv.EulerByFactorization(tc.n)
			resFou := dv.EulerByFourier(tc.n)

			if resDef != tc.expected {
				t.Errorf("По определению: для n=%d ожидалось %d, а получили %d", tc.n, tc.expected, resDef)
			}
			if resFac != tc.expected {
				t.Errorf("По факторизации: для n=%d ожидалось %d, а получили %d", tc.n, tc.expected, resFac)
			}
			if resFou != tc.expected {
				t.Errorf("Через Фурье: для n=%d ожидалось %d, а получили %d", tc.n, tc.expected, resFou)
			}
		})
	}
}

func TestRootsOfUnity(t *testing.T) {
	dv := DataVerifier{}

	tests := []int{1, 2, 3, 4, 8}

	for _, n := range tests {
		roots := dv.RootsOfUnity(n)

		if len(roots) != n {
			t.Errorf("Для n=%d ожидалось %d корней, получено %d", n, n, len(roots))
		}

		// (z^n == 1)
		for _, root := range roots {

			powered := cmplx.Pow(root, complex(float64(n), 0))

			if cmplx.Abs(powered-1) > 1e-9 {
				t.Errorf("Ошибка: корень %v в степени %d дал %v (а должен был дать 1)", root, n, powered)
			}
		}
	}
}

func TestPrimitiveRootsOfUnity(t *testing.T) {
	dv := DataVerifier{}

	tests := []int{3, 4, 5, 6, 8}

	for _, n := range tests {
		primitives := dv.PrimitiveRootsOfUnity(n)

		// кол-во корней = результат функции Эйлера
		expectedCount := dv.EulerByDef(n)
		if len(primitives) != expectedCount {
			t.Errorf("Для n=%d ожидалось %d первообразных корней, получено %d", n, expectedCount, len(primitives))
		}

		// проверка свойств корня
		for _, root := range primitives {

			// в степени n даёт 1
			poweredN := cmplx.Pow(root, complex(float64(n), 0))
			if cmplx.Abs(poweredN-1) > 1e-9 {
				t.Errorf("Ошибка: корень %v в степени %d не равен 1", root, n)
			}

			// от 1 до n-1 он не даёт 1
			for k := 1; k < n; k++ {
				poweredK := cmplx.Pow(root, complex(float64(k), 0))

				if cmplx.Abs(poweredK-1) < 1e-9 {
					t.Errorf("Корень %v фальшивый! Он дал 1 в степени %d (а должен был только в %d)", root, k, n)
				}
			}
		}
	}
}

func TestVandermondeMatrices(t *testing.T) {
	dv := DataVerifier{}
	n := 4 // 4x4

	// первообразный корень
	roots := dv.PrimitiveRootsOfUnity(n)
	omega := roots[0]

	// матрица и обратная матрица
	V, Vinv := dv.VandermondeMatrices(n, omega)

	// перемножение матриц
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			var sum complex128 = 0
			for k := 0; k < n; k++ {
				sum += V[i][k] * Vinv[k][j]
			}

			// проверяемм элементы нна диагонали (1)
			if i == j {
				if cmplx.Abs(sum-1) > 1e-9 {
					t.Errorf("На диагонали в ячейке [%d][%d] ожидалась 1, получено %v", i, j, sum)
				}
			} else {
				// вне дигонали (0)
				if cmplx.Abs(sum) > 1e-9 {
					t.Errorf("В ячейке [%d][%d] ожидался 0, получено %v", i, j, sum)
				}
			}
		}
	}
}

func TestFourierTransform(t *testing.T) {
	dv := DataVerifier{}

	tests := []struct {
		name string
		vec  []complex128
	}{
		{"Обычные числа", []complex128{1, 2, 3, 4}},
		{"Массив нулей", []complex128{0, 0, 0, 0}},
		{"Числа с минусом", []complex128{10, -5, 3, 7}},
		{"Комплексные числа", []complex128{complex(1, 1), complex(2, -2)}},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {

			originalVector := tc.vec
			n := len(originalVector)

			// первообразный корень для размера массива
			primitives := dv.PrimitiveRootsOfUnity(n)
			omega := primitives[0] // берем первый

			// получаем матрицы
			V, Vinv := dv.VandermondeMatrices(n, omega)

			// прямое преобразование
			transformed := dv.FurierTransform(V, originalVector)

			// обратное
			restoredVector := dv.FurierTransform(Vinv, transformed)

			for i := 0; i < n; i++ {
				// cравниваем исходный массив с восстановленным
				if cmplx.Abs(originalVector[i]-restoredVector[i]) > 1e-9 {
					t.Errorf("Ошибка Фурье! Вектор не восстановился.\nБыло: %v\nСтало: %v",
						originalVector[i], restoredVector[i])
				}
			}
		})
	}
}
