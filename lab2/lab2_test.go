package lab2

import (
	"reflect"
	"testing"
)

// нахождения первообразных корней
func TestFindAllPrimitiveRoots(t *testing.T) {
	dv := DataVerifier{}

	tests := []struct {
		name          string
		n             int
		expectRoots   bool // сущеествование корней
		expectedCount int  // колличество корней phi(phi(n))
	}{
		{"простое 5", 5, true, 2},
		{"9 (степень 3)", 9, true, 2},
		{"8 (корней нет)", 8, false, 0},
		{"12 (корней нет)", 12, false, 0},
		{"простое 17", 17, true, 8},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			roots := dv.FindAllPrimitiveRoots(tc.n)

			// должны быть корни?
			if tc.expectRoots && len(roots) == 0 {
				t.Errorf("Для n=%d ожидались корни, но функция вернула пустой список", tc.n)
			}
			if !tc.expectRoots && len(roots) > 0 {
				t.Errorf("Для n=%d корней быть не должно, но они найдены: %v", tc.n, roots)
			}

			// сколько их?
			if len(roots) != tc.expectedCount {
				t.Errorf("Для n=%d ожидалось %d корней, получено %d", tc.n, tc.expectedCount, len(roots))
			}

			// первообразный ли он?
			if len(roots) > 0 {
				phi := dv.euler(tc.n)
				for _, root := range roots {
					// степень phi -> 1 по модулю n
					if power(root, phi, tc.n) != 1 {
						t.Errorf("Корень %d не первообразный: в степени %d не дал 1", root, phi)
					}
					// в меньшей степени не 1
					factors := primeFactors(phi)
					for p := range factors {
						if power(root, phi/p, tc.n) == 1 {
							t.Errorf("Корень %d зациклился слишком рано (в степени %d)", root, phi/p)
						}
					}
				}
			}
		})
	}
}

// матрица Вандермонда
func TestVandermondeMatrixMod(t *testing.T) {
	dv := DataVerifier{}
	n := 17
	m := 4
	omega := 13
	omegaInv := modInverse(omega, n)

	t.Run("прямая матрица на обратную", func(t *testing.T) {
		// матрицы
		V := dv.VandermondeMatrixMod(m, n, omega)
		VinvUnscaled := dv.VandermondeMatrixMod(m, n, omegaInv)

		mInv := modInverse(m, n)

		// V на Vinv
		for i := 0; i < m; i++ {
			for j := 0; j < m; j++ {
				sum := 0
				for k := 0; k < m; k++ {
					// перемножаем и делим на m (по модулю n)
					term := (V[i][k] * VinvUnscaled[k][j]) % n
					term = (term * mInv) % n
					sum = (sum + term) % n
				}

				// единичная матрица
				if i == j && sum != 1 {
					t.Errorf("На диагонали ожидалась 1, получено %d", sum)
				}
				if i != j && sum != 0 {
					t.Errorf("Вне диагонали ожидался 0, получено %d", sum)
				}
			}
		}
	})
}

// медленный ДПФ
func TestSlowFurierTrans(t *testing.T) {
	dv := DataVerifier{}

	n := 17
	m := 4
	omega := 13
	omegaInv := modInverse(omega, n)
	mInv := modInverse(m, n)

	tests := []struct {
		name string
		arr  []int
	}{
		{"Обычные числа", []int{1, 2, 3, 4}},
		{"Нули", []int{0, 0, 0, 0}},
		{"Одинаковые числа", []int{5, 5, 5, 5}},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			original := tc.arr

			V := dv.VandermondeMatrixMod(m, n, omega)
			VinvUnscaled := dv.VandermondeMatrixMod(m, n, omegaInv)

			// прямое преобразование
			transformed := dv.SlowFurierTrans(V, original, n)

			// обратное
			restoredUnscaled := dv.SlowFurierTrans(VinvUnscaled, transformed, n)

			restored := make([]int, m)
			for i := 0; i < m; i++ {
				restored[i] = (restoredUnscaled[i] * mInv) % n
			}

			// сравниваем
			if !reflect.DeepEqual(original, restored) {
				t.Errorf("Ошибка! Ожидалось: %v, получено: %v", original, restored)
			}
		})
	}
}

// быстрое преобразование Фурьье
func TestFastAFurierTrans(t *testing.T) {
	dv := DataVerifier{}

	n := 17
	omega := 13 // корень степени 4 по модулю 17

	tests := []struct {
		name string
		arr  []int
	}{
		{"Обычные числа", []int{1, 2, 3, 4}},
		{"Нули", []int{0, 0, 0, 0}},
		{"Одинаковые числа", []int{5, 5, 5, 5}},
		{"Возрастающие числа", []int{10, 11, 12, 13}},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			original := tc.arr

			// прямое быстрое преобразование
			transformed := dv.FastFurierTrans(original, n, omega)

			// обратное
			restored := dv.InvFastFurierTrans(transformed, n, omega)

			// cравниваем
			if !reflect.DeepEqual(original, restored) {
				t.Errorf("Ошибка! Ожидалось: %v, получено: %v", original, restored)
			}
		})
	}
}
