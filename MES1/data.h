#pragma once
#include<iostream>
#include<cmath>
using namespace std;


class Wezel {
private:
	int x; //wspolrzedna wezla
	int war_brz; //warunek brzegowy 0 - brak, 1 - strumien, 2-konwekcja
	double alfa;
	double q;
	double t_sr;
public:
	Wezel() {

	}
	Wezel(int a, int b, double alfa, double q, double t_sr) {
		x = a;
		war_brz = b;
		this->alfa = alfa;
		this->q = q;
		this->t_sr = t_sr;
	}
	void set_x(int a) {
		x = a;
	}
	int get_x() {
		return x;
	}
	int get_war_brz() {
		return war_brz;
	}

	double get_alfa() {
		return alfa;
	}
	double get_q() {
		return q;
	}
	double get_t_sr() {
		return t_sr;
	}
};

class Element {
private:
	int id_el;
	int id_1;
	int id_2;
	double L;
	double s;
	double k;
	double** H;
	double* P;
	double C;
public:
	Element(int id_el, int id_1, int id_2, double l, double s, double k);
	~Element();
	int get_id_el() {
		return id_el;
	}
	int get_id_1() {
		return id_1;
	}
	int get_id_2() {
		return id_2;
	}
	double get_L() {
		return L;
	}
	double get_s() {
		return s;
	}
	double get_k() {
		return k;
	}
	double** get_H() {
		return H;
	}
	double* get_P() {
		return P;
	}

	void create(double alfa_S, double q_S, double alfa_tsr_S);
	void wysw_lokal_H();
	void wysw_lokal_P();
};

class Siatka {
private:
	int l_el;
	int l_wez;
	double **global_H;
	double *global_P;
	double *temp;

	double **HP;

	Wezel** tab_wezl;
	Element** tab_el;
public:
	Siatka(int l_wez);
	~Siatka();
	int get_l_wez() {
		return l_wez;
	}
	int get_l_el() {
		return l_el;
	}
	double** get_global_H() {
		return global_H;
	}
	double* get_global_P() {
		return global_P;
	}
	double* get_temp() {
		return temp;
	}

	void generuj_siatke(double l, double s, double k, double alfa, double q, double t_sr);
	void wysw_lokalne_macierze_elementow();
	void generuj_H_P();
	void wysw_global_H_P();
	void wysw_HP();

	void gauss();

	void wysw_temp();
};
