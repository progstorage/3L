#ifdef _WIN32
	#include "stdafx.h"
#endif

#include "config.h"
#include "source/3L.cpp"
#include "source/Menu.cpp"


int main(void) {
	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len = ME.x.size();
	double d;	//"диаметр" фигуры - максимальное расстояние между точками, 0.04*d - расстояние до точек уровней 

	system("python rmdir.py");	//удаляем старые файлы 
	fs::create_directory("sorted_points");	//создаем папку, куда сохранятся файлы с группами точек
	fs::create_directory("levels");		//создаем папку, куда сохранятся файлы с уровнями для каждой группы
	fs::create_directory("coefs");		//создаем папку, куда сохранятся файлы с коэффициентами многочлена для каждой группы

	switch (ME.execution_method) {
	case 1:
		old_sort(ME.x, ME.y, &d);
		break;
	case 2:
		sort_method_2(ME.x, ME.y, &d);
		break;
	default:
		cerr << "Nonexistent menu method" << endl;
		exit(0);
	}

	vector<double> x;
	vector<double> y;
	for (int i = 0; i < points_clusters_array[0].size(); i++) {
		x.push_back(points_clusters_array[0][i].first);
		y.push_back(points_clusters_array[0][i].second);
	}
	cubic_spline Spline;
	Spline.build_spline(x, y, 11);
	//Spline.write_poinst(Spline, x[0], x[x.size()-1]);
	double res = Spline.f(0);
	cout << "res = " << res;
	cin.get(); cin.get();

	setup();

	fill_levels_multithread(d);
	system("python plot/plot_3L.py");	// отрисовка

	return 0;
}

