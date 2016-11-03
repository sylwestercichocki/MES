//#include<iostream>
#include"data.h"
//using namespace std;


int main() {
	int l_w=0;
	double L=0;
	double s=0;
	double k=0;
	double alfa = 0;
	double q = 0;
	double t_sr = 0;

	cout << "Podaj ilosc wezlow" << endl;
	cin >> l_w;
	cout << "Podaj dlugosc preta" << endl;
	cin >> L;
	cout << "Podaj pole przekroju s" << endl;
	cin >> s;
	cout << "Podaj k" << endl;
	cin >> k;
	cout << "Podaj alfa" << endl;
	cin >> alfa;
	cout << "Podaj q" << endl;
	cin >> q;
	cout << "Podaj temperature otoczenia" << endl;
	cin >> t_sr;


	Siatka siatka(l_w);
	siatka.generuj_siatke(L, s, k, alfa, q, t_sr);
	cout << endl;
	siatka.wysw_lokalne_macierze_elementow();
	cout << endl;
	siatka.generuj_H_P();
	siatka.wysw_global_H_P();
	siatka.wysw_HP();

	siatka.gauss();
	siatka.wysw_HP();
	siatka.wysw_temp();



	system("PAUSE");
	return 0;
}
