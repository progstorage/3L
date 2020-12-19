#ifndef MENU
#define MENU

#include <vector>

class Menu {
public:
	int execution_method; //1 - ñòàðûé ìåòîä ñîåäèíåíèÿ òî÷åê, 2 - íîâûé ìåòîä ñ ðàçäåëåíèåì òî÷åê íà ãðóïïû
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
	void select_execution_method();
};

#endif