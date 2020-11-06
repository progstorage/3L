#include "./source/3L.h"

int main(void) {

	string path = "tests/in/";
	Menu ME(path); //Меню выбора тесового файла с координатами точек

	int len;
	double d;	//"диаметр" фигуры - максимальное расстояние между точками, 0.04*d - расстояние до точек уровней
	vector_pairs sorted_points;
	char plot[50];
	
	//while (ME.to_draw) {
	len = ME.x.size();
	sorted_points = sort(ME.x, ME.y, &d);

	Matrix M(len * 3, vector<double>(6));
	vector_pairs l1(len);
	vector_pairs l3(len);

	levels(sorted_points, d, l1, l3);
	fill_3L_Matrix_2nd_power(M, sorted_points, l1, l3); 
	solve_system(M, d);
	//Sleep(2000);
	ME.test_file = "python plot/plot_3L.py " + ME.test_file;
		
	strcpy_s(plot, ME.test_file.c_str());
	system(plot);		// отрисовка
	
						
	//ME.start();

	// clear variables
	/*len = 0;
	d = 0;
	sorted_points.clear();
	M.~vector();
	l1.~vector();
	l3.~vector();*/
	//}
	
	return 0;
}