#pragma once
#include<vector>
#include<iostream>
					//Oh God, the comments are in GridPack.cpp Im going asleep
struct Node {
	double x, y;	//Koordy nod�w ||| x to szeroko��, a y to wysoko��
	short int bc;	//Is border node or not
	Node(double x1, double y1, short int bc1);

	void print();
};

class Element {		
public:
	int id[4];		//Id Nod�w w rogach elementu (Nody poi�dzy kt�rymi jest element)
	double h[4][4];
	double jacobian[2][2];
	double jacobianInversed[2][2];
	double jacobDet;
	double jacobDetInversed;
	double(*derivNEta)[4];
	double(*derivNKsi)[4];
	int nIntegrationPoints;
	double hbc[4][4];
	double p[4];
	double k1;		//Wspolczynnik do liczenia h (dla materialu)
	double alpha;	//Wspolczynnik konwekcji do liczenia hbc
	double tSurrounding; //Temperatrua otoczenia

	Element(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, int points, double(*derivNEta)[4], double(*derivNKsi)[4], double alpha1, double hbc1[4][4], double tSurrounding1, double p1[4]);
	//~Element();

	static Element* create9(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, double alpha1, double tSurrounding);

	static Element* create4(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, double alpha1, double tSurrounding);

	//void initialization(double k, std::vector<Node> nodes);

	void printids();

	int* returnIdsArrayPointer();

	void printNEta();

	void printNKsi();

	void printH();

	void printHbc();

	void printP();

};


class Grid {

public:
	double H;	//Wysoko��
	double B;	//Szeroko��
	int nH;		//Ilo�� nod�w w wysoko�ci
	int nB;		//Ilo�� nod�w w szeroko�ci
	int nN;		//Ilo�� wszystkich nod�w
	int nE;		//Ilo�� wszystkich element�w
	std::vector<std::vector<double>> hGlobal; //One H to rule them all, One H to find them, One H to bring them all and in the element bind them. (Potezne H calej siatki)
	std::vector<double> pGlobal; //P as above for p

	double deltaH;	//Wysoko�� jednego elementu
	double deltaB;	//Szeroko�� jednego elementu

	std::vector<Node> nodes;		//Fajniej dynamiczna tablica nod�w
	std::vector<Element> elements;	//Fajniej dynamiczna tablica element�w

	Grid(double H1, double B1, int nH1, int nB1, int n = 4, double k = 30.0, double alpha1 = 25.0, double tSurrounding = 1200.0);

	void printAllElements();

	void printEntireGrid();

	void printHGlobal();

	void printPGlobal();


};