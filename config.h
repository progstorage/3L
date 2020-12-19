#pragma once
#include <omp.h>
/**
 * ����� ������������ OMP (���������� ��� ������������� ����������):
 * ������ -> �������� -> C/C++ -> ���� -> ��������� Open MP -> �� (���� ������������ Visual Studio)
 * ���� ���, �� ������� �������� ���� ����������: icc -openmp <>.cpp
 */

// �����:
// - 0: �������
// - 1: ��������
// - 2: �������� � �������� ���������
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