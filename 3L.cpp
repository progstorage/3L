#include "3L.h"


/*-----------------------------------------------------------------------------------------------------------------
---------------------------------------------------NEW ALG---------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------*/
vector<double> matr_max(const Matrix& M) {
	// возвращает вектор максимальных значений каждого столбца матрицы 
	int cols = M.cols_num();
	int rows = M.rows_num();
	vector<double> res;
	double tmpmax;
	for (int i = 0; i < cols; i++) {
		tmpmax = M[0][i];
		for (int j = 0; j < rows; j++) {
			if (M[j][i] > tmpmax) {
				tmpmax = M[j][i];
			}
		}
		res.push_back(tmpmax);
	}
	return res;
}

int vec_min_ind(vector<double> v) {
	// индекс наименьшего элемента вектора
	double min = v[0];
	double min_ind = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[i] < min) {
			min = v[i];
			min_ind = i;
		}
	}
	return min_ind;
}

double vec_max(vector<double> v) {
	// максимальный элемент вектора
	double max = v[0];
	for (int i = 0; i < v.size(); i++) {
		if (v[i] > max) {
			max = v[i];
		}
	}
	return max;
}

double vec_sum(vector<double> v, int deg) {
	// сумма значений вектора, возведенных в степень deg
	double sum = 0;
	for (int i = 0; i < v.size(); i++) {
		sum += pow(v[i], deg);
	}
	return sum;
}

vector<double> vec_mult(vector<double> v1, vector<double> v2) {
	// поэлементное умножение векторов 
	vector<double> res;
	for (int i = 0; i < v1.size(); i++) {
		res.push_back(v1[i] * v2[i]);
	}
	return res;
}

double mean(vector<double> x) {
	// среднее значение вектора
	double res = 0;
	for (int i = 0; i < x.size(); i++) {
		res += x[i];
	}
	res /= x.size();
	return res;
}


void quickSort(vector<double>& arr, int left, int right, vector<int>& ind) {
	int i = left, j = right;
	double tmp;
	int tmpind;
	double pivot = arr[(left + right) / 2];
	while (i <= j) {
		while (arr[i] < pivot)
			i++;
		while (arr[j] > pivot)
			j--;
		if (i <= j) {
			tmp = arr[i];
			arr[i] = arr[j];
			arr[j] = tmp;

			tmpind = ind[i];
			ind[i] = ind[j];
			ind[j] = tmpind;

			i++;
			j--;
		}
	};
	if (left < j)
		quickSort(arr, left, j, ind);
	if (i < right)
		quickSort(arr, i, right, ind);
}


void mls(double* x0, double* y0, double* d, vector<double> x, vector<double> y) {
	// moving least squares
	// задана точка(x0, y0) и набор её соседей(x(i), y(i))
	vector<double> xx = x;
	vector<double> yy = y;
	xx.insert(xx.begin(), *x0);
	yy.insert(yy.begin(), *y0);
	// центр тяжести
	double xc = mean(xx);
	double yc = mean(yy);

	int k = 16;
	vector<double> p;
	for (int i = 0; i < k; i++) {
		p.push_back(i * 2 * 3.14159 / k);
	}
	//находим угол, при повороте на который отклонение точек от новой оси x минимально
	vector<double> tmpx, tmpy, sin_p, cos_p, tmp;
	for (int i = 0; i < xx.size(); i++) {
		tmpx.push_back(xx[i] - xc);
		tmpy.push_back(yy[i] - yc);
	}
	for (int i = 0; i < p.size(); i++) {
		sin_p.push_back(sin(p[i]));
		cos_p.push_back(cos(p[i]));
	}

	Matrix M(tmpx.size(), sin_p.size());

	for (int i = 0; i < p.size(); i++) {
		for (int j = 0; j < xx.size(); j++) {
			M[j][i] = abs(tmpx[j] * sin_p[i] + tmpy[j] * cos_p[i]);
		}
	}
	int min_ind = vec_min_ind(matr_max(M));

	for (int i = 0; i < xx.size(); i++) {
		tmp.push_back(tmpx[i] * cos(p[min_ind]) - tmpy[i] * sin(p[min_ind]));
		yy[i] = tmpx[i] * sin(p[min_ind]) + tmpy[i] * cos(p[min_ind]);
	}
	xx = tmp;
	// в повёрнутых координатах аппроксимируем параболой
	double a11, a12, a13, a23, a33, b1, b2, b3;
	a11 = xx.size();
	a12 = vec_sum(xx, 1);
	a13 = vec_sum(xx, 2);
	a23 = vec_sum(xx, 3);
	a33 = vec_sum(xx, 4);

	b1 = vec_sum(yy, 1);
	b2 = vec_sum(vec_mult(xx, yy), 1);
	b3 = vec_sum(vec_mult(vec_mult(xx, xx), yy), 1);

	Matrix A(3, 3);
	A[0][0] = a11; A[0][1] = a12; A[0][2] = a13;
	A[1][0] = a12; A[1][1] = a13; A[1][2] = a23;
	A[2][0] = a13; A[2][1] = a23; A[2][2] = a33;

	Matrix B(3, 1);
	B[0][0] = b1;
	B[1][0] = b2;
	B[2][0] = b3;

	Matrix C(3, 1);
	C = A.inverse()*B;
	// сглаженные координаты первой точки и ширина коридора
	Matrix D(xx.size(), 3);

	for (int i = 0; i < xx.size(); i++) {
		D[i][0] = 1;
		D[i][1] = xx[i];
		D[i][2] = xx[i] * xx[i];
	}
	Matrix Ans = D * C;
	vector<double> ans;
	for (int i = 0; i < xx.size(); i++) {
		ans.push_back(abs(Ans[i][0] - yy[i]));
	}
	*d = vec_max(ans);

	double xs = xx[0];
	Matrix Tmp(1, 3);
	Tmp[0][0] = 1; Tmp[0][1] = xs; Tmp[0][2] = xs * xs;
	double ys = (Tmp*C)[0][0];
	// обратное преобразование координат
	double tmpd = xs * cos(p[min_ind]) + ys * sin(p[min_ind]) + xc;
	ys = -xs * sin(p[min_ind]) + ys * cos(p[min_ind]) + yc;
	xs = tmpd;

	*x0 = xs;
	*y0 = ys;
}


vector<double> pre3L(vector<double>& x, vector<double>& y) {
	// предварительное сглаживание данных
	int n = 8;
	vector<double> xs, ys, d;
	for (int i = 0; i < x.size(); i++) {
		xs.push_back(0);
	}
	ys = xs; d = xs;

	for (int i = 0; i < x.size(); i++) {
		vector<double> r2, tmpx, tmpy;
		vector<int> tmpind, tmpind_2;
		for (int j = 0; j < x.size(); j++) {
			r2.push_back((x[j] - x[i])*(x[j] - x[i]) + (y[j] - y[i])*(y[j] - y[i]));
			tmpind.push_back(j);
		}
		quickSort(r2, 0, r2.size() - 1, tmpind);
		for (int k = 0; k < tmpind.size(); k++) {
			if (r2[k] != 0) {
				tmpind_2.push_back(tmpind[k]);
			}
			if (tmpind_2.size() == n) {
				break;
			}
		}

		xs[i] = x[i]; ys[i] = y[i];
		for (int e = 0; e < tmpind_2.size(); e++) {
			tmpx.push_back(x[tmpind_2[e]]);
			tmpy.push_back(y[tmpind_2[e]]);
		}
		mls(&xs[i], &ys[i], &d[i], tmpx, tmpy);
	}
	// отбрасываем точки, для которых "коридор" втрое шире среднего
	double mean_d = mean(d);
	vector<double> xs_new, ys_new, d_new;
	for (int i = 0; i < d.size(); i++) {
		if (d[i] <= 3 * mean_d) {
			xs_new.push_back(xs[i]);
			ys_new.push_back(ys[i]);
			d_new.push_back(d[i]);
		}
	}
	xs = xs_new; ys = ys_new; d = d_new;
	// пары очень близких точек заменяем на одиночные точки
	double r2, flag = 0;
	vector<double> resx, resy;
	for (int i = 0; i < xs.size(); i++) {
		for (int j = 0; j < x.size() - 1; j++) {
			r2 = (xs[j + 1] - xs[i])*(xs[j + 1] - xs[i]) + (ys[j + 1] - ys[i])*(ys[j + 1] - ys[i]);
			if (r2 > 9 * mean_d * mean_d) {
				resx.push_back(xs[i]); resy.push_back(ys[i]);
			}
		}
	}
	x = xs; y = ys;

	return d;
}



/*-----------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------*/






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



void test_spline_norm_levels(vector_pairs points, double d, vector_pairs& l1, vector_pairs& l3, int num) {
	// TEST
	// Построение уровней с помощью нормалей к сплайнам, построенным по точкам points
	string type = ".txt";
	string dir = "levels/";

	int len = points.size();

	vector<double> x;
	vector<double> y;
	for (int i = 0; i < points_clusters_array[0].size(); i++) {
		x.push_back(points[i].first);
		y.push_back(points[i].second);
	}
	cubic_spline Spline;
	Spline.build_spline(x, y, len);
	double nx[3], ny[3], dy, n;

	for (int i = 0; i < len; i++) {
		nx[0] = x[i];
		nx[1] = x[i] + d;
		nx[2] = nx[1] - nx[0];

		if (i != 0) {
			dy = Spline.diff_spline(x[i]);
		}
		else {
			// в первой точке производная считается неправильно, поэтому отступаем от нее на малую величину
			dy = Spline.diff_spline(x[i] + 0.000001);
		}

		ny[0] = Spline.f(nx[0]);
		ny[1] = Spline.f(nx[1]) - d / dy;
		ny[2] = ny[1] - ny[0];

		n = sqrt(nx[2] * nx[2] + ny[2] * ny[2]);
		nx[2] /= n;
		ny[2] /= n;

		if (dy > 0) {
			l1[i].first = x[i] + d * nx[2];
			l1[i].second = y[i] + d * ny[2];
			l3[i].first = x[i] - d * nx[2];
			l3[i].second = y[i] - d * ny[2];
		}
		else {
			l1[i].first = x[i] - d * nx[2];
			l1[i].second = y[i] - d * ny[2];
			l3[i].first = x[i] + d * nx[2];
			l3[i].second = y[i] + d * ny[2];
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

	Matrix a;

	#if TEST_MODE == 1
		clock_t start_1, end_1, start_2, end_2;
	
		// многопоточное умножение матрицы на вектор
		//auto start_1 = std::chrono::system_clock::now();
		start_1 = clock();
		a = ((M.transpose() * M).inverse() * M.transpose()) * b.transpose();
		end_1 = clock() - start_1;

		//auto end_1 = std::chrono::duration<double>(std::chrono::system_clock::now() - start_1);
			
		// однопоточное
		//auto start_2 = std::chrono::system_clock::now();
		start_2 = clock();
		Matrix _a = mult_Matrix(mult_Matrix(mult_Matrix(M.transpose(), M).inverse(), M.transpose()), b.transpose());
		end_2 = clock() - start_2;
		//auto end_2 = std::chrono::duration<double>(std::chrono::system_clock::now() - start_2);
		
		cout << "\nSingle-threaded multiplication execution cpu time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution cpu time:\t" << end_1 / (double)CLOCKS_PER_SEC;
	
		if (_a == a)
			cout << "\nBoth methods gave the same result!" << endl;
		else
			cout << "\nDifferent results, something went wrong!" << endl;
	#elif TEST_MODE == 2
		// big-size matrix testing
		clock_t start_1, end_1, start_2, end_2;

		// многопоточное умножение матрицы на вектор
		//auto start_1 = std::chrono::system_clock::now();
		// a = mult_Matrix_multithread(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		start_1 = clock();
		a = generate_random_Matrix(size_1, size_2) * generate_random_Matrix(size_2, size_3);
		end_1 = clock() - start_1;
		//auto end_1 = std::chrono::duration<double>(std::chrono::system_clock::now() - start_1);


		// однопоточное
		//auto start_2 = std::chrono::system_clock::now();
		start_2 = clock();
		Matrix _a = mult_Matrix(generate_random_Matrix(size_1, size_2), generate_random_Matrix(size_2, size_3));
		end_2 = clock() - start_2;
		//auto end_2 = std::chrono::duration<double>(std::chrono::system_clock::now() - start_2);

		cout << "\nSingle-threaded multiplication execution cpu time:\t" << end_2 / (double)CLOCKS_PER_SEC
			 << "\nMulti-threaded multiplication execution cpu time:\t" << end_1 / (double)CLOCKS_PER_SEC;
		cout << endl;
	#else
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
	Matrix M;

	#if TEST_MODE
		// TEST_MODE 1 or 2
		clock_t start_1, end_1, start_2, end_2;
		//auto start_1 = std::chrono::system_clock::now();
		start_1 = clock();
		// многопоточное создание уровней
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
		//auto end_1 = std::chrono::duration<double>(std::chrono::system_clock::now() - start_1);
		end_1 = clock() - start_1;
		
		// однопоточное
		//auto start_2 = std::chrono::system_clock::now();
		start_2 = clock();
		for (int j = 0; j < points_clusters_count; j++) {
			l1.resize(points_clusters_array[j].size());
			l3.resize(points_clusters_array[j].size());
			levels(points_clusters_array[j], d, l1, l3, j + 1);

			M.resize(points_clusters_array[j].size() * 3, 6);
			M.fill_3L_Matrix_2nd_power(points_clusters_array[j], l1, l3);
			solve_system(M, d, j + 1);

			l1.clear();
			l3.clear();
		}
		//auto end_2 = std::chrono::duration<double> (std::chrono::system_clock::now() - start_2);
		end_2 = clock() - start_2;
		cout << "\nSingle-threaded level creation cpu time:\t" << end_2 / (double)CLOCKS_PER_SEC
			<< "\nMulti-threaded level creation cpu time:\t" << end_1 / (double)CLOCKS_PER_SEC;

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