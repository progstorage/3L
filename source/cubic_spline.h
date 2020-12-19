#ifndef CUBIC_SPLINE
#define CUBIC_SPLINE

#include <vector>

class cubic_spline {
private:
	// Ñòðóêòóðà, îïèñûâàþùàÿ ñïëàéí íà êàæäîì ñåãìåíòå ñåòêè
	struct spline_tuple {
		double a, b, c, d, x;
	};
	spline_tuple *splines; // Ñïëàéí
	std::size_t n; // Êîëè÷åñòâî óçëîâ ñåòêè
	void free_mem(); // Îñâîáîæäåíèå ïàìÿòè
public:
	cubic_spline(); //êîíñòðóêòîð
	~cubic_spline(); //äåñòðóêòîð
	// Ïîñòðîåíèå ñïëàéíà
	// x - óçëû ñåòêè, äîëæíû áûòü óïîðÿäî÷åíû ïî âîçðàñòàíèþ, êðàòíûå óçëû çàïðåùåíû
	// y - çíà÷åíèÿ ôóíêöèè â óçëàõ ñåòêè
	// n - êîëè÷åñòâî óçëîâ ñåòêè
	void build_spline(vector<double>&, vector<double>&, std::size_t);
	void write_poinst(cubic_spline, double, double, int); // Çàïèñûâàåò çíà÷åíèÿ ñïëàéíà â òî÷êàõ ïðîìåæóòêà [xmin-3, xmax+3]
	double f(double) const;	// Âû÷èñëåíèå çíà÷åíèÿ èíòåðïîëèðîâàííîé ôóíêöèè â ïðîèçâîëüíîé òî÷êå
};

#endif