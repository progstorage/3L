#include "cubic_spline.h"

cubic_spline::cubic_spline() : splines(NULL) { }
cubic_spline::~cubic_spline() { free_mem(); }

void cubic_spline::build_spline(vector<double>& x, vector<double>& y, std::size_t n) {
	free_mem();
	this->n = n;
	// Инициализация массива сплайнов
	splines = new spline_tuple[n];
	for (std::size_t i = 0; i < n; ++i) {
		splines[i].x = x[i];
		splines[i].a = y[i];
	}
	splines[0].c = 0.;
	// Решение СЛАУ относительно коэффициентов сплайнов c[i] методом прогонки для трехдиагональных матриц
	// Вычисление прогоночных коэффициентов - прямой ход метода прогонки
	double *alpha = new double[n - 1];
	double *beta = new double[n - 1];
	double A, B, C, F, h_i, h_i1, z;
	alpha[0] = beta[0] = 0.;
	for (std::size_t i = 1; i < n - 1; ++i) {
		h_i = x[i] - x[i - 1], h_i1 = x[i + 1] - x[i];
		A = h_i;
		C = 2 * (h_i + h_i1);
		B = h_i1;
		F = 6 * ((y[i + 1] - y[i]) / h_i1 - (y[i] - y[i - 1]) / h_i);
		z = (A * alpha[i - 1] + C);
		alpha[i] = -B / z;
		beta[i] = (F - A * beta[i - 1]) / z;
	}
	splines[n - 1].c = (F - A * beta[n - 2]) / (C + A * alpha[n - 2]);
	// Нахождение решения - обратный ход метода прогонки
	for (std::size_t i = n - 2; i > 0; --i)
		splines[i].c = alpha[i] * splines[i + 1].c + beta[i];
	// Освобождение памяти, занимаемой прогоночными коэффициентами
	delete[] beta;
	delete[] alpha;
	// По известным коэффициентам c[i] находим значения b[i] и d[i]
	for (std::size_t i = n - 1; i > 0; --i) {
		double h_i = x[i] - x[i - 1];
		splines[i].d = (splines[i].c - splines[i - 1].c) / h_i;
		splines[i].b = h_i * (2. * splines[i].c + splines[i - 1].c) /
			6 + (y[i] - y[i - 1]) / h_i;
	}
}

double cubic_spline::f(double x) const {// Возвращает значение интерполированной функции в точке x
	if (!splines)
		return std::numeric_limits<double>::quiet_NaN(); // Если сплайны ещё не построены - возвращаем NaN
	spline_tuple *s;
	if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл - тов массива
		s = splines + 1;
	else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл - том массива
		s = splines + n - 1;
	else { // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл - та массива
		std::size_t i = 0, j = n - 1;
		while (i + 1 < j) {
			std::size_t k = i + (j - i) / 2;
			if (x <= splines[k].x)
				j = k;
			else
				i = k;
		}
		s = splines + j;
	}
	double dx = (x - s->x);
	return s->a + (s->b + (s->c / 2 + s->d * dx / 6.) * dx) * dx; // Вычисляем значение сплайна в заданной точке по схеме Горнера
}

void cubic_spline::free_mem()
{
	delete[] splines;
	splines = NULL;
}

void cubic_spline::write_poinst(cubic_spline Spline, double min, double max, int numpoints = 100) {
	vector<double> x;
	vector<double> y;
	int i = 0;
	double p = min;
	double h = (max - min) / numpoints;
	while (p <= max) {
		x.push_back(p);
		y.push_back(Spline.f(p));
		p += h;
	}

	ofstream outfile("out.txt", ios::out | ios::trunc); //хз че за ошибка, нажимаешь продолжить и все работает
	for (int i = 0; i <= numpoints; i++) {
		outfile << fixed << x[i] << "," << y[i] << endl;
	}
	outfile.close();
}


void cubic_spline::write_poinst(cubic_spline Spline, double min, double max, int numpoints = 100) {
	vector<double> x;
	vector<double> y;
	int i = 0;
	double p = min;
	double h = (max - min) / numpoints;
	while (p <= max) {
		x.push_back(p);
		y.push_back(Spline.f(p));
		p += h;
	}

	ofstream outfile("out.txt", ios::out | ios::trunc); //хз че за ошибка, нажимаешь продолжить и все работает
	for (int i = 0; i <= numpoints; i++) {
		outfile << fixed << x[i] << "," << y[i] << endl;
	}
	outfile.close();
}


double cubic_spline::diff_spline(double x) const {// Возвращает значение интерполированной функции в точке x
	if (!splines)
		return std::numeric_limits<double>::quiet_NaN(); // Если сплайны ещё не построены - возвращаем NaN
	spline_tuple *s;
	if (x <= splines[0].x) // Если x меньше точки сетки x[0] - пользуемся первым эл - тов массива
		s = splines + 1;
	else if (x >= splines[n - 1].x) // Если x больше точки сетки x[n - 1] - пользуемся последним эл - том массива
		s = splines + n - 1;
	else { // Иначе x лежит между граничными точками сетки - производим бинарный поиск нужного эл - та массива
		std::size_t i = 0, j = n - 1;
		while (i + 1 < j) {
			std::size_t k = i + (j - i) / 2;
			if (x <= splines[k].x)
				j = k;
			else
				i = k;
		}
		s = splines + j;
	}
	double dx = (x - s->x);
	double res = s->b + s->c + dx + s->d * dx * dx / 2.0;
	if (res) {
		return res;
	}
	else {
		return 0.0000000000000000001;
	}
	//return s->a + (s->b + (s->c / 2 + s->d * dx / 6.) * dx) * dx; // Вычисляем значение сплайна в заданной точке по схеме Горнера
}