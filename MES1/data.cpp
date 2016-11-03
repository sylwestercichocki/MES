#include "data.h"

Siatka::Siatka(int l_wez)
{
	this->l_wez = l_wez;
	l_el = l_wez - 1;

	global_H = new  double*[l_wez];
	for (int i = 0; i < l_wez; ++i) {
		global_H[i] = new double[l_wez];
	}

	HP = new double*[l_wez];

	for (int i = 0; i < l_wez; ++i) {
		HP[i] = new double[l_wez + 1];
	}

	global_P = new double[l_wez];
	temp = new double[l_wez];

	for (int i = 0; i < l_wez; ++i) {
		for (int j = 0; j < l_wez; ++j) {
			global_H[i][j] = 0;
		}
		global_P[i] = 0;
		temp[i] = 0;
	}

	for (int i = 0; i < l_wez; ++i) {
		for (int j = 0; j < l_wez+1; ++j) {
			HP[i][j] = 0;
		}
	}
}

Siatka::~Siatka()
{
	for (int i = 0; i < l_wez; ++i) {
		delete[] global_H[i];
	}
	delete[] global_H;
	delete[] global_P;
	delete[] temp;
	delete[] tab_wezl;
	delete[] tab_el;
}

void Siatka::generuj_siatke(double L, double s, double k, double alfa, double q, double t_sr)
{
	double l = L / static_cast<double>(l_el);
	tab_wezl = new Wezel*[l_wez];
	tab_el = new Element*[l_el];
	for (int i = 0; i < l_wez; ++i) {
		if (i == 0) tab_wezl[i] = new Wezel(i, 1, alfa, q, t_sr);
		else if (i == l_wez-1) tab_wezl[i] = new Wezel(i, 2, alfa, q, t_sr);
		else tab_wezl[i] = new Wezel(i,0, 0 ,0 ,0);
	}
	for (int i = 0; i < l_el; ++i) {
		tab_el[i] = new Element(i, i, i + 1, l, s, k);
	}


	for (int i = 0; i < l_el; ++i) {
		if (tab_wezl[tab_el[i]->get_id_1()]->get_war_brz() == 1) {		//dopisuje odpowiednie warunki
			tab_el[i]->create(0, q*tab_el[i]->get_s(), 0);			//brzegowe do macierzy lokalnych
		}
		else if (tab_wezl[tab_el[i]->get_id_2()]->get_war_brz() == 2) {
				tab_el[i]->create(alfa*tab_el[i]->get_s(), 0, -alfa*t_sr*tab_el[i]->get_s());
		}
		else {
			tab_el[i]->create(0, 0, 0);
		}
		
	}

}

void Siatka::wysw_lokalne_macierze_elementow()
{
	for (int i = 0; i < l_el; ++i) {
		tab_el[i]->wysw_lokal_H();
		tab_el[i]->wysw_lokal_P();
	}
}

void Siatka::generuj_H_P()
{
	
	for (int i = 0; i < l_el; ++i) {
		
		double **h = tab_el[i]->get_H();
		int ii = tab_el[i]->get_id_el();
		global_H[ii][ii] += h[0][0];
		global_H[ii][ii+1] += h[0][1];
		global_H[ii+1][ii] += h[1][0];
		global_H[ii+1][ii+1] += h[1][1];
		
	}

	for (int i = 0; i < l_el; ++i) {
		double *p = tab_el[i]->get_P();
		global_P[tab_el[i]->get_id_el()] += p[0];
		global_P[tab_el[i]->get_id_el()+1] += p[1];
	}


	for (int i = 0; i < l_wez; ++i) {
		for (int j = 0; j < l_wez; ++j) {
			HP[i][j] = global_H[i][j];
		}
	}

	for (int i = 0; i < l_wez; ++i) {
		HP[i][l_wez] = global_P[i];
	}
}

void Siatka::wysw_global_H_P()
{
	cout << "Globalna macierz H" << endl;
	for (int i = 0; i < l_wez; ++i) {
		for (int j = 0; j < l_wez; ++j) {
			cout << global_H[i][j] << " ";
		}
		cout << endl;
	}

	cout << "Globalna macierz P" << endl;
	for (int i = 0; i < l_wez; ++i) {
		cout << global_P[i] << endl;
	}

	cout << endl;
}

void Siatka::wysw_HP()
{
	for (int i = 0; i < l_wez; ++i) {
		for (int j = 0; j < l_wez + 1; ++j) {
			cout << HP[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
}

void Siatka::gauss()
{

	for (int i = 0; i<l_wez; i++) {
		// Search for maximum in this column
		double maxEl = abs(HP[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<l_wez; k++) {
			if (abs(HP[k][i]) > maxEl) {
				maxEl = abs(HP[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<l_wez + 1; k++) {
			double tmp = HP[maxRow][k];
			HP[maxRow][k] = HP[i][k];
			HP[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<l_wez; k++) {
			double c = -HP[k][i] / HP[i][i];
			for (int j = i; j<l_wez + 1; j++) {
				if (i == j) {
					HP[k][j] = 0;
				}
				else {
					HP[k][j] += c * HP[i][j];
				}
			}
		}
	}
	// Solve equation Ax=b for an upper triangular matrix A
	for (int i = l_wez - 1; i >= 0; i--) {
		temp[i] = -(HP[i][l_wez] / HP[i][i]);
		for (int k = i - 1; k >= 0; k--) {
			HP[k][l_wez] += HP[k][i] * temp[i];
		}
	}

}

void Siatka::wysw_temp()
{
	cout << "Koncowe temperatury" << endl;
	for (int i = 0; i < l_wez; ++i) {
		cout << "T" << i << ": " << temp[i] << endl;
	}
}

Element::Element(int id_el, int id_1, int id_2, double l, double s, double k)
{
	this->id_el = id_el;
	this->id_1 = id_1;
	this->id_2 = id_2;
	this->L = l;
	this->s = s;
	this->k = k;
	C = (s*k) / L;
	H = new double*[2];
	for (int i = 0; i < 2; ++i) {
		H[i] = new double[2];
	}
	P = new double[2];

	H[0][0] = 0;
	H[0][1] = 0;
	H[1][0] = 0;
	H[1][1] = 0;

	P[0] = 0;
	P[1] = 0;
}

Element::~Element()
{
	for (int i = 0; i < 2; ++i) {
		delete[] H[i];
	}
	delete[] H;
	delete[] P;
}

void Element::create(double alfa_S, double q_S, double alfa_tsr_S) //generuje macierze lokalne
{																   //dla kazdego elementu
	H[0][0] = C;
	H[0][1] = -C;
	H[1][0] = -C;
	H[1][1] = C + alfa_S;

	P[0] = q_S;
	P[1] = alfa_tsr_S;
}

void Element::wysw_lokal_H()
{
	cout << "Macierz lokalna H elementu nr: " << id_el << endl;
	cout << H[0][0] << " " << H[0][1] << endl;
	cout << H[1][0] << " " << H[1][1] << endl;
}

void Element::wysw_lokal_P()
{
	cout << "Macierz lokalna P elementu nr: " << id_el << endl;
	cout << P[0] << endl << P[1] << endl;
}
