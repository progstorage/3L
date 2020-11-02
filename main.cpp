#include "./source/3L.h"

int main(void) {
	string path = "tests/in/";
	Menu ME(path); //���� ������ �������� ����� � ������������ �����

	int len = ME.x.size();
	double d;	//"�������" ������ - ������������ ���������� ����� �������, 0.04*d - ���������� �� ����� ������� 
	vector_pairs sorted_point = sort(ME.x, ME.y, &d);

	vector_pairs l1(len);
	vector_pairs l3(len);
	levels(sorted_point, d, l1, l3);

	Matrix M(len * 3, vector<double>(6));
	fill_3L_Matrix_2nd_power(M, sorted_point, l1, l3);

	solve_system(M, d);

	while (ME.to_draw) {
		ME.test_file = "python plot/plot_3L.py " + ME.test_file;
		char plot[50];
		strcpy_s(plot, ME.test_file.c_str());
		system(plot);		// ���������
		ME.start();
	}

	sorted_point.clear();

	return 0;
}