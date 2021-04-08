#ifndef MENU
#define MENU

#include <vector>

class Menu {
public:
	int matrix_power; // степень матрицы
	string dir;
	vector<string> files;
	int n;
	vector<double> x;
	vector<double> y;
	string test_file;
	bool to_draw;

	Menu(string);
	void get_files(string);
	void print_menu();
	void select_file();
	void split(string, double*, double*);
	void read_file();
	void start();
	void select_matrix_power();
};

#endif