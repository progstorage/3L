#include "3L.h"

double diff(pair<double, double> pk, pair<double, double> pk_1) {
	//производная в точке pk = (x[k], y[k]) 
	//pk_1 = (x[k+1], y[k+1])
	return ((pk_1.second - pk.second) / (pk_1.first - pk.first));
}


int old_sort(const vector<double>& x, const vector<double>& y, double* d) {
	// функция сортирует исходные координаты в смысле евклидовой метрики
	// и возвращает вектор пар

	int len = x.size();

	vector_pairs sorted_points(len);
	sorted_points[0].first = x[0];
	sorted_points[0].second = y[0];
	vector<int> flags(len);

	for (int i = 0; i < len; i++) {
		flags[i] = 0;
	}
	double dist; //расстояние до первой точки
	double min_dist = sqrt((x[1] - x[0])*(x[1] - x[0]) + (y[1] - y[0])*(y[1] - y[0]));
	double max_dist = min_dist;
	int min_dist_ind = 0;
	double dydx;
	int i = 0;
	for (int k = 0; k < len - 1; k++) {
		if (flags[i] == 0) {
			for (int j = 0; j < len; j++) {
				if ((j != i) && (flags[j] == 0)) {
					dist = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
					if (dist >= max_dist) { max_dist = dist; }
				}
			}
			sorted_points[k + 1].first = x[min_dist_ind];
			sorted_points[k + 1].second = y[min_dist_ind];
			flags[i] = 1;
			i = min_dist_ind;
			min_dist = 100;
		}
	}
	*d = max_dist * 0.04;

	points_clusters_array.push_back(sorted_points);//отсортированная группа точек добавляется в массив групп
	points_clusters_count++;

	ofstream outfile("sorted_points/sorted_points.txt", ios::out | ios::trunc);
	for (int i = 0; i < len; i++) {
		outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
	}

	return 0;
}


int sort_method_2(const vector<double>& x, const vector<double>& y, double* d) {
	// функция сортирует исходные координаты в смысле евклидовой метрики
	// если расстояние до следующей точки "слишком большое" (в данном случае - больше расстояния до первой точки),
	// то все точки разбиваются на отдельные группы, в которых расстояния от каждой точки до следующей не превышает этого значения
	// => формируется глобальный вектор отсортированных векторов пар, а также каждая группа точек сохраняется в файл {число}.txt

	*d = 0;
	int len = x.size();

	vector_pairs sorted_points(len);
	sorted_points[0].first = x[0];
	sorted_points[0].second = y[0];
	vector<int> flags(len);

	for (int i = 0; i < len; i++) {
		flags[i] = 0;
	}
	double dist, dist_to_fitrs;
	double min_dist = sqrt((x[1] - x[0])*(x[1] - x[0]) + (y[1] - y[0])*(y[1] - y[0]));
	double max_dist = min_dist;
	int min_dist_ind = 0;
	int i = 0;
	int retflag;
	string type = ".txt";
	string dir = "sorted_points/";
	for (int k = 0; k < len - 1; k++) {
		if (flags[i] == 0) {
			flags[i] = 1;
			for (int j = 0; j < len; j++) {
				if ((j != i) && (flags[j] == 0)) {
					retflag = 1;//если нашлась точка, в которой еще не были
					dist = sqrt((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
					if (dist <= min_dist) {
						min_dist = dist;
						min_dist_ind = j;
					}
					if (dist >= max_dist) { max_dist = dist; }
				}
			}
			dist_to_fitrs = sqrt((x[min_dist_ind] - x[0])*(x[min_dist_ind] - x[0]) + (y[min_dist_ind] - y[0])*(y[min_dist_ind] - y[0]));
			if (min_dist <= dist_to_fitrs) {
				sorted_points[k + 1].first = x[min_dist_ind]; 
				sorted_points[k + 1].second = y[min_dist_ind];
				i = min_dist_ind;
				min_dist = 100;
				retflag = 0;
			}
			else {//если расстояние до следующей точки больше, чем расстояние до первой точки (вообще говоря, если больше диаметра группы,
				  //но как определить диаметр группы пока не понятно
				vector_pairs tmp_sort_points(k + 1);
				for (int j = 0; j < k + 1; j++) {
					tmp_sort_points[j] = sorted_points[j];
				}
				points_clusters_array.push_back(tmp_sort_points);//отсортированная группа точек добавляется в массив групп
				points_clusters_count++;
				auto name = dir + to_string(points_clusters_count) + type;
				ofstream outfile(name, ios::out | ios::trunc);
				for (int i = 0; i < k + 1; i++) {
					outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
				}

				flags[i] = 1;
				flags[min_dist_ind] = 1;
				vector<double> new_x;//оставшиеся точки формируют новые массивы координат х и у
				vector<double> new_y;
				for (int i = 0; i < len; i++) {
					if (flags[i] == 0) {
						new_x.push_back(x[i]);
						new_y.push_back(y[i]);
					}
				}

				if (new_x.size() != 0) { sort_method_2(new_x, new_y, d); }//рекурсивный вызов функции для оставшихся точек
				else { return 0; }
			}
		}
	}
	if (retflag != 1) {//если перебраны все точки
		points_clusters_array.push_back(sorted_points);
		points_clusters_count++;
		auto name = dir + to_string(points_clusters_count) + type;
		ofstream outfile(name, ios::out | ios::trunc);
		for (int i = 0; i < sorted_points.size(); i++) {
			outfile << fixed << sorted_points[i].first << "," << sorted_points[i].second << endl;
		}
		return 0;
	}
	*d = max_dist * 0.04;
	return 0;
}


void levels(vector_pairs points, double d, vector_pairs& l1, vector_pairs& l3, int num) {
	// функция принимает на вход пары точек
	// и записывает полученные уровни в txt файл
	string type = ".txt";
	string dir = "levels/";

	int len = points.size();

	double tmpx, tmpy, n;
	for (int i = 1; i < len - 1; i++) {
		tmpx = points[i + 1].second - points[i - 1].second;
		tmpy = -(points[i + 1].first - points[i - 1].first);
		n = sqrt((tmpx*tmpx + tmpy*tmpy));
		tmpx = tmpx / n;
		tmpy = tmpy / n;
		l1[i].first = points[i].first + d * tmpx;
		l1[i].second = points[i].second + d * tmpy;
		l3[i].first = points[i].first - d * tmpx;
		l3[i].second = points[i].second - d * tmpy;
		if (i == 1) {
			l1[0].first = points[0].first + d * tmpx;
			l1[0].second = points[0].second + d * tmpy;
			l3[0].first = points[0].first - d * tmpx;
			l3[0].second = points[0].second - d * tmpy;
		}
		if (i == len - 2) {
			l1[len - 1].first = points[len - 1].first + d * tmpx;
			l1[len - 1].second = points[len - 1].second + d * tmpy;
			l3[len - 1].first = points[len - 1].first - d * tmpx;
			l3[len - 1].second = points[len - 1].second - d * tmpy;
		}
	}

	auto name = dir + to_string(num) + type;
	ofstream outfile(name, ios::out | ios::trunc);
	for (int i = 0; i < len; i++) {
		outfile << fixed << l1[i].first << "," << l1[i].second << "," << l3[i].first << "," << l3[i].second << endl;
	}
	outfile.close();
}


Matrix generate_random_Matrix(const int n, const int m) {
	// генерирует случайную матрицу размера n x m
	Matrix res(n, m);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			res[i][j] = (rand() % 10000 - 5000) / 100;
	return res;
}


void solve_system(Matrix& M, double d, int num) {
	// Решает систему уравнений - находит коэффициенты многочлена и записывает их в файл "coef.txt"
	string type = ".txt";
	string dir = "coefs/";
	double len = M.rows_num();
	Matrix b(1, len);
	for (int i = 0; i < len / 3; i++) { b[0][i] = d; }
	for (int i = len / 3; i < 2 * len / 3; i++) { b[0][i] = 0; }
	for (int i = 2 * len / 3; i < len; i++) { b[0][i] = -d; }

	Matrix a(0, 0); // возможно ошибка (не делал конструктор по умолчанию и не перегружал оператор =)

	#if TEST_MODE == 1
		clock_t start_1, end_1, start_2, end_2;
	
		// многопоточное умножение матрицы на вектор
		start_1 = clock();
		// a = mult_Matrix_multithread(mult_Matrix_multithread(inverse_Matrix(mult_Matrix_multithread(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));		
		a = ((M.transpose() * M).inverse() * M.transpose()) * b.transpose();
		end_1 = clock() - start_1;
			
		// однопоточное
		start_2 = clock();
		Matrix _a = mult_Matrix(mult_Matrix(inverse_Matrix(mult_Matrix(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));
		end_2 = clock() - start_2;
		
		cout << "\nSingle-threaded multiplication execution time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution time:\t" << end_1 / (double)CLOCKS_PER_SEC;
	
		if (_a == a)
			cout << "\nBoth methods gave the same result!" << endl;
		else
			cout << "\nDifferent results, something went wrong!" << endl;
	#elif TEST_MODE == 2
		// big-size matrix testing
		clock_t start_1, end_1, start_2, end_2;

		// многопоточное умножение матрицы на вектор
		start_1 = clock();
		// a = mult_Matrix_multithread(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		a = generate_random_Matrix(size_1, size_2) * generate_random_Matrix(size_2, size_3);
		end_1 = clock() - start_1;

		// однопоточное
		start_2 = clock();
		Matrix _a = mult_Matrix(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		end_2 = clock() - start_2;

		cout << "\nSingle-threaded multiplication execution time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution time:\t" << end_1 / (double)CLOCKS_PER_SEC;
		cout << endl;
	#else
		// a = mult_Matrix_multithread(mult_Matrix_multithread(inverse_Matrix(mult_Matrix_multithread(transpose_Matrix(M), M)), transpose_Matrix(M)), transpose_Matrix(b));
		a = ((M.transpose() * M).inverse() * M.transpose()) * b.transpose();
	#endif

	auto name = dir + to_string(num) + type;
	ofstream outfile(name, ios::out | ios::trunc);
	for (int i = 0; i < 6; i++) {
		outfile << fixed << a[i][0] << endl;
	}
	outfile.close();
}


void fill_levels_multithread(double d) {
	int i;
	vector_pairs l1, l3;
	Matrix M(0, 0);

	#if TEST_MODE
		clock_t start_1, end_1, start_2, end_2;
		start_1 = clock();

		// многопоточное создание уровней
		#pragma omp parallel for private(i) shared(points_clusters_array, l1, l3)
			for (i = 0; i < points_clusters_count; i++) {
				l1.resize(points_clusters_array[i].size());
				l3.resize(points_clusters_array[i].size());

				levels(points_clusters_array[i], d, l1, l3, i + 1);

				M.resize(points_clusters_array[i].size() * 3, 6);
				fill_3L_Matrix_2nd_power(M, points_clusters_array[i], l1, l3);
				solve_system(M, d, i + 1);

				l1.clear();
				l3.clear();
			}
		end_1 = clock() - start_1;

		// однопоточное
		start_2 = clock();
		for (int j = 0; j < points_clusters_count; j++) {
			l1.resize(points_clusters_array[j].size());
			l3.resize(points_clusters_array[j].size());
			levels(points_clusters_array[j], d, l1, l3, j + 1);

			M.resize(points_clusters_array[j].size() * 3, 6);
			fill_3L_Matrix_2nd_power(M, points_clusters_array[j], l1, l3);
			solve_system(M, d, j + 1);

			l1.clear();
			l3.clear();
		}
		end_2 = clock() - start_2;

		cout << "\nSingle-threaded level creation time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded level creation time:\t" << end_1 / (double)CLOCKS_PER_SEC;

	#else
		#pragma omp parallel for private(i) shared(points_clusters_array, l1, l3)
			for (i = 0; i < points_clusters_count; i++) {
				l1.resize(points_clusters_array[i].size());
				l3.resize(points_clusters_array[i].size());

				levels(points_clusters_array[i], d, l1, l3, i + 1);

				M.resize(points_clusters_array[i].size() * 3, 6);
				M.fill_3L_Matrix_2nd_power(points_clusters_array[i], l1, l3);
				solve_system(M, d, i + 1);

				l1.clear();
				l3.clear();
			}
	#endif
}