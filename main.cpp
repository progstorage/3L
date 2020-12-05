//#include "stdafx.h"
#include "source/3L.cpp"

/**
 * Чтобы использовать OMP (библиотека для многопоточных вычислений):
 * Проект -> Свойства -> C/C++ -> Язык -> Поддержка Open MP -> Да (если используется Visual Studio)
 * если нет, то следует добавить ключ компиляции: icc -openmp <>.cpp	
 */

int main(void) {
	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len = ME.x.size();
	double d;	//"диаметр" фигуры - максимальное расстояние между точками, 0.04*d - расстояние до точек уровней 

	system("python rmdir.py");//удаляем старые файлы 
	fs::create_directory("sorted_points");//создаем папку, куда сохранятся файлы с группами точек
	fs::create_directory("levels");//создаем папку, куда сохранятся файлы с уровнями для каждой группы
	fs::create_directory("coefs");//создаем папку, куда сохранятся файлы с коэффициентами многочлена для каждой группы

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

	const int threadsNum_mult = 6;
	omp_set_dynamic(0);
	omp_set_nested(1);
	omp_set_num_threads(threadsNum_mult);

	fill_levels_multithread(d);

	system("python plot/plot_3L.py");	// отрисовка

	return 0;
}

