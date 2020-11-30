// Project.cpp: определяет точку входа для консольного приложения.
//

//#include "stdafx.h"
#include "./source/3L.cpp"

int main(void) {
	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len = ME.x.size();
	double d;	//"диаметр" фигуры - максимальное расстояние между точками, 0.04*d - расстояние до точек уровней 

	system("python rmdir.py");//удаляем старые файлы 
	fs::create_directory("sorted_points");//создаем папку, куда сохранятся файлы с группами точек
	fs::create_directory("levels");//создаем папку, куда сохранятся файлы с уровнями для каждой группы
	fs::create_directory("coefs");//создаем папку, куда сохранятся файлы с коэффициентами многочлена для каждой группы

	if (ME.execution_method == 1) {
		old_sort(ME.x, ME.y, &d);
	}
	else {
		sort_method_2(ME.x, ME.y, &d);
	}
	
	for (int i = 0; i < points_clusters_count; i++) {
		vector_pairs l1(points_clusters_array[i].size());
		vector_pairs l3(points_clusters_array[i].size());
		levels(points_clusters_array[i], d, l1, l3, i+1);

		Matrix M(points_clusters_array[i].size() * 3, vector<double>(6));
		fill_3L_Matrix_2nd_power(M, points_clusters_array[i], l1, l3);
		solve_system(M, d, i + 1);

		l1.clear();
		l3.clear();
	}
	
	system("python plot/plot_3L.py");	// отрисовка

	return 0;
}

