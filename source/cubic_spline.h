#include <vector>

class cubic_spline
{
private:
	// Ñòðóêòóðà, îïèñûâàþùàÿ ñïëàéí íà êàæäîì ñåãìåíòå ñåòêè
	struct spline_tuple
	{
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
	void build_spline(vector<double>& x, vector<double>& y, std::size_t n);
	void write_poinst(cubic_spline Spline, double min, double max, int numpoints); // Çàïèñûâàåò çíà÷åíèÿ ñïëàéíà â òî÷êàõ ïðîìåæóòêà [xmin-3, xmax+3]
	double f(double x) const;	// Âû÷èñëåíèå çíà÷åíèÿ èíòåðïîëèðîâàííîé ôóíêöèè â ïðîèçâîëüíîé òî÷êå
};