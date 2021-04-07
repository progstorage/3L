#ifdef _WIN32
//#include "stdafx.h"
#endif

#include "config.h"
#include "source/3L.cpp"
#include "source/Menu.cpp"

//#include "source/3L.h"
//#include "source/Matrix.cpp"

int main(void) {
	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len = ME.x.size();
	//double d;	//"диаметр" фигуры - максимальное расстояние между точками, 0.04*d - расстояние до точек уровней 

	system("python rmdir.py");	//удаляем старые файлы 
	fs::create_directory("sorted_points");	//создаем папку, куда сохранятся файлы с группами точек
	fs::create_directory("levels");		//создаем папку, куда сохранятся файлы с уровнями для каждой группы
	fs::create_directory("coefs");		//создаем папку, куда сохранятся файлы с коэффициентами многочлена для каждой группы

	double d; 
	vector<Matrix> res_3L;
	switch (ME.execution_method) {
	case 1:
		d = pre3L(ME.x, ME.y);
		res_3L = alg3L(ME.x, ME.y);
		writeToFile(res_3L);
		//old_sort(ME.x, ME.y, &d);
		break;
	case 2:
		//vector<double> d = pre3L(ME.x, ME.y);
		d = pre3L(ME.x, ME.y);
		res_3L = alg3L(ME.x, ME.y);
		writeToFile(res_3L);
		// sort_method_2(ME.x, ME.y, &d);
		break;
	default:
		cerr << "Nonexistent menu method" << endl;
		exit(0);
	}

	/*for (int i = 0; i < ME.x.size(); i++) {
		cout << ME.x[i] << " ";
	}*/

	// vector<double> x;
	// vector<double> y;
	// for (int i = 0; i < points_clusters_array[0].size(); i++) {
	// 	x.push_back(points_clusters_array[0][i].first);
	// 	y.push_back(points_clusters_array[0][i].second);
	// }
	// cubic_spline Spline;
	// Spline.build_spline(x, y, 11);
	// //Spline.write_poinst(Spline, x[0], x[x.size()-1]);
	// double res = Spline.f(0);
	//cout << "res = " << res;
	//cin.get(); cin.get();

	setup();
	double d1 = 0.008;
	int num = 1;
	vector_pairs l1(ME.x.size()), l3(ME.x.size()), points(ME.x.size());
	for (int i = 0; i < ME.x.size(); i++) {
		points[i].first = ME.x[i]; 
		points[i].second = ME.y[i];
	}

	levels(points, d, l1, l3, num);


	/*double len_l1 = l1.size();
	Matrix M(len_l1 *3, 6);
	M.fill_3L_Matrix_2nd_power(points, l1, l3);

	solve_system(M, d1);*/


	//fill_levels_multithread(d1);
	cout << "\nPlotting...";
	system("python plot/plot_3L.py");	// отрисовка

	return 0;
}

