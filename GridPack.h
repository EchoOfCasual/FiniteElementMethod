#pragma once
#include<vector>
#include<iostream>
#include "Elements_2D.h"

struct Node {
	double x, y;	//Koordy nodów ||| x to szerokoœæ, a y to wysokoœæ
	Node(double x1, double y1);

	void printCoords();
};

class Element {
public:
	int id[4];		//Id Nodów w rogach elementu (Nody poiêdzy którymi jest element)
	double h[4][4];
	double jacobian[2][2];
	double jacobianInv[2][2];
	double jacobDet;
	double jacobDetInv;
	double derivNEta[1][1];
	double derivNKsi[1][1];

	Element(int id1 = 0, int id2 = 0 , int  id3 = 0, int id4 = 0);


	void printids();

	int* returnIdsArrayPointer();
};


class Grid {

public:
	double H;	//Wysokoœæ
	double B;	//Szerokoœæ
	int nH;		//Iloœæ nodów w wysokoœci
	int nB;		//Iloœæ nodów w szerokoœci
	int nN;		//Iloœæ wszystkich nodów
	int nE;		//Iloœæ wszystkich elementów

	double deltaH;	//Wysokoœæ jednego elementu
	double deltaB;	//Szerokoœæ jednego elementu

	std::vector<Node> nodes;		//Fajniej dynamiczna tablica nodów
	std::vector<Element> elements;	//Fajniej dynamiczna tablica elementów

	Grid(double H1, double B1, int nH1, int nB1, int n = 4, double k = 30.0);

	void printAllElements();

	void printEntireGrid();

	void calculate(int n=4, double k=30.0);

	void jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element4_2D elem);
	void jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element9_2D elem);
};