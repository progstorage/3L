// 3L_main.cpp: определяет точку входа для консольного приложения.
//

#include "stdafx.h"

#include "config.h"
#include "source/3L.cpp"
#include "source/Menu.cpp"

int main(void) {
	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len = ME.x.size();

	system("python rmdir.py");	//удаляем старые файлы 
	fs::create_directory("sorted_points");	//создаем папку, куда сохранятся файлы с группами точек
	fs::create_directory("levels");		//создаем папку, куда сохранятся файлы с уровнями для каждой группы
	fs::create_directory("coefs");		//создаем папку, куда сохранятся файлы с коэффициентами многочлена для каждой группы

	double d;
	vector<Matrix> res_3L;
	d = pre3L(ME.x, ME.y);
	res_3L = alg3L(ME.x, ME.y);
	writeToFile(res_3L);

	vector_pairs points;
	for (Matrix mp : res_3L) {
		for (int i = 0; i < mp.rows_num(); i++) {
			points.push_back({ mp[i][0], mp[i][1] });
		}
	}
	setup();

	vector_pairs l1(points.size()), l3(points.size());
	levels(points, d, l1, l3); // построение уровней

	double len_l1 = l1.size();
	Matrix M;
	int max_col = 0;
	for (int i = 1; i <= ME.matrix_power + 1; i++) {
		max_col += i;
	}
	M.resize(len_l1 * 3, max_col);
	M.fill_3L_Matrix_power(points, l1, l3, ME.matrix_power);
		
	solve_system(M, d); // решение СЛАУ с помощью псевдообратной матрицы

	//fill_levels_multithread(d1);
	cout << "\nPlotting...";
	system("python plot/plot_3L.py");	// отрисовка
	cout << "\nEnd";
	cin.get(); cin.get();
	return 0;
}

