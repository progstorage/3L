#pragma once
#include <omp.h>
/**
 * Чтобы использовать OMP (библиотека для многопоточных вычислений):
 * Проект -> Свойства -> C/C++ -> Язык -> Поддержка Open MP -> Да (если используется Visual Studio)
 * если нет, то следует добавить ключ компиляции: icc -openmp <>.cpp
 */

// режим:
// - 0: обычный
// - 1: тестовый
// - 2: тестовый с большими матрицами
#define TEST_MODE 0

#if TEST_MODE == 2
	// M(size_1 x size_2), N(size_2 x size_3), res(size_1, size_3)
	const int size_1 = 100;
	const int size_2 = 20;
	const int size_3 = 400;
#endif

// multi-thread conf
const int threadsNum_mult = 6;

void setup() {
	omp_set_dynamic(0);
	omp_set_nested(1);
	omp_set_num_threads(threadsNum_mult);
}