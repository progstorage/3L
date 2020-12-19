#include "Menu.h"


Menu::Menu(string path_to_dir) {
	dir = path_to_dir;
	get_files(path_to_dir);
	start();
	to_draw = true;
}


void Menu::get_files(string path_to_dir) {
	for (const auto & f : fs::directory_iterator(path_to_dir)) {
		files.push_back(f.path().string());
	}
}


void Menu::print_menu() {
	cout << "Directory: " << dir << endl;
	cout << "File" << "\t\t\t" << "N" << endl;
	
	int i;
	for (i = 0; i < files.size(); i++) {
		cout << files[i] << "\t" << i + 1 << endl;
	}

	cout << "Exit" << "\t\t\t" << i + 1 << endl;
}


void Menu::select_file() {
	cout << endl << "Enter file number:\t";
	cin >> n;
	if (n <= files.size()) {
		test_file = files[n - 1];
	}
	else if (n == files.size() + 1) {
		to_draw = false;
		exit(0);
	}
	else {
		cout << endl << "Error! Try again!";
		select_file();
	}
	cout << "\n-------------------------\n\n";
}

void Menu::select_execution_method() {
	cout << endl << "Enter 1 - old method, 2 - new method:\t";
	cin >> execution_method;
}

void Menu::start() {
	print_menu();
	select_file();
	read_file();
	select_execution_method();
}


void Menu::split(string line, double* x, double* y) {
	auto pos = line.find(";");
	if (pos != string::npos) {
		*x = stod(line.substr(0, pos));
		*y = stod(line.substr(pos + 1));
	}
}


void Menu::read_file() {
	string line;
	double tmpx;
	double tmpy;
	ifstream in(test_file);		//Чтение данных из файла test_file
	if (in.is_open()) {
		while (getline(in, line)) {
			split(line, &tmpx, &tmpy);
			x.push_back(tmpx);
			y.push_back(tmpy);
		}
	}
	in.close();
}