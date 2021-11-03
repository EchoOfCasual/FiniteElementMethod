#include <iostream>
#include "Gauss.h"
#include "Elements_2D.h"
#include <vector>
#include "GridPack.h"


void jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element9_2D elem, Grid grid) {
	//The very same imp belllowww, not tested yet tho, soooooooo, good luck with that future me
	jakobian[0][0] = 0;
	jakobian[0][1] = 0;
	jakobian[1][0] = 0;
	jakobian[1][1] = 0;

	//double xHardcode[4] = { 0, 0.025, 0.025, 0 };
	//double yHardcode[4] = { 0, 0, 0.025, 0.025 };

	for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
		jakobian[1][0] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobian[1][1] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
	}

	/*for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * xHardcode[i];
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * yHardcode[i];
		jakobian[1][0] += elem.derivNEta[pCJ][i] * xHardcode[i];
		jakobian[1][1] += elem.derivNEta[pCJ][i] * yHardcode[i];
	}*/


	jakobianInv[0][0] = 0;
	jakobianInv[0][1] = 0;
	jakobianInv[1][0] = 0;
	jakobianInv[1][1] = 0;

	for (int i = 0; i < 4; i++) {
		jakobianInv[1][1] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobianInv[0][1] -= elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
		jakobianInv[1][0] -= elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobianInv[0][0] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
	}

};

void jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element4_2D elem, Grid grid) {
	jakobian[0][0] = 0;
	jakobian[0][1] = 0;
	jakobian[1][0] = 0;
	jakobian[1][1] = 0;

	//double xHardcode[4] = { 0, 0.025, 0.025, 0 };
	//double yHardcode[4] = { 0, 0, 0.025, 0.025 };

	for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
		jakobian[1][0] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobian[1][1] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
	}

	/*for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * xHardcode[i];
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * yHardcode[i];
		jakobian[1][0] += elem.derivNEta[pCJ][i] * xHardcode[i];
		jakobian[1][1] += elem.derivNEta[pCJ][i] * yHardcode[i];
	}*/


	jakobianInv[0][0] = 0;
	jakobianInv[0][1] = 0;
	jakobianInv[1][0] = 0;
	jakobianInv[1][1] = 0;

	for (int i = 0; i < 4; i++) {
		jakobianInv[1][1] += elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobianInv[0][1] -= elem.derivNKsi[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
		jakobianInv[1][0] -= elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].x;
		jakobianInv[0][0] += elem.derivNEta[pCJ][i] * grid.nodes[grid.elements[elemI].id[i]].y;
	}


};

class Sth {
public:
	double jakob[2][2];
	double jakobInversed[2][2];
	int nPC;
	Grid grid = Grid(0.025, 0.025, 2, 2);

	Sth(Grid g = Grid(0.025, 0.025, 2, 2), int n = 4) {
		grid = g;
		nPC = n;

		for (int i = 0; i < grid.nE; i++) {
			for (int j = 0; j < nPC; j++) {
				if (nPC == 4) {
					jakobian(i, j, jakob, jakobInversed, Element4_2D(), grid);
				}
				else if (nPC == 9) {
					jakobian(i, j, jakob, jakobInversed, Element9_2D(), grid);
				}
				else {
					std::cout << "Nope.";
				}

				double det = jakob[0][0] * jakob[1][1] - jakob[1][0] * jakob[0][1];
				double detInv = 1.0 / det;

				std::cout << "\ndet: "<< det << " ||| detInv: " << detInv << "\n============================================\n";
				for (int d = 0; d < 2; d++) {
					for (int e = 0; e < 2; e++) {
						std::cout << jakob[d][e] << " \t\t ";
					}
					std::cout << "\n";
				}
				std::cout << "--------------------------------------------\n";
				double mutiplied[2][2];

				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						mutiplied[i][j] = jakobInversed[i][j] * detInv;
						std::cout << mutiplied[i][j] << "\t\t";
					}
					std::cout << "\n";
				}
				std::cout << "\n\n";





			}
		}


	}


};


/*class Element9_2D {
public:
	double derivNEta[9][4];
	double derivNKsi[9][4];

	Element9_2D() {
		for (int i = 0; i < 9; i++) {

			double eta, ksi;

			if (i == 0 || i == 3 || i == 6) {
				ksi = -sqrt(3.0/5.0);
			}
			else if (i == 1 || i == 4 || i == 7) {
				ksi = 0;
			}
			else {
				ksi = sqrt(3.0/ 5.0);
			}

			if (i == 0 || i == 1 || i == 2) {
				eta = -sqrt(3.0 / 5.0);
			}
			else if (i == 3 || i == 4 || i == 5) {
				eta = 0;
			}
			else {
				eta = sqrt(3.0 / 5.0);
			}

			derivNEta[i][0] = -0.25 * (1 - ksi);
			derivNEta[i][1] = -0.25 * (1 + ksi);
			derivNEta[i][2] = 0.25 * (1 + ksi);
			derivNEta[i][3] = 0.25 * (1 - ksi);

			derivNKsi[i][0] = -0.25 * (1 - eta);
			derivNKsi[i][1] = 0.25 * (1 - eta);
			derivNKsi[i][2] = 0.25 * (1 + eta);
			derivNKsi[i][3] = -0.25 * (1 + eta);



		}
	}

	void printNEta() {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNEta[i][j] << "\t";
				if (i == 1 || i == 4 || i == 7) {
					std::cout << "\t";
				}
			}
			std::cout << "\n";
		}
	}

	void printNKsi() {
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNKsi[i][j] << "\t";
				
				if (i > 2 && i < 6) {
					std::cout << "\t";
				}
				
			}
			std::cout << "\n";
		}
	}
};

class Element4_2D {
public:
	double derivNEta[4][4];
	double derivNKsi[4][4];

	Element4_2D() {
		for (int i = 0; i < 4; i++) {
			double eta, ksi;
			if (i == 0 || i == 1) {
				eta = -1 / sqrt(3);
			}
			else {
				eta = 1 / sqrt(3);
			}

			if (i == 0 || i == 3) {
				ksi = -1 / sqrt(3);
			}
			else {
				ksi = 1 / sqrt(3);
			}

			derivNEta[i][0] = -0.25 * (1 - ksi);
			derivNEta[i][1] = -0.25 * (1 + ksi);
			derivNEta[i][2] = 0.25 * (1 + ksi);
			derivNEta[i][3] = 0.25 * (1 - ksi);

			derivNKsi[i][0] = -0.25 * (1 - eta);
			derivNKsi[i][1] = 0.25 * (1 - eta);
			derivNKsi[i][2] = 0.25 * (1 + eta);
			derivNKsi[i][3] = -0.25 * (1 + eta);



		}
	}

	void printNEta() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNEta[i][j] << "\t";
			}
			std::cout << "\n";
		}
	}

	void printNKsi() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNKsi[i][j] << "\t";
			}
			std::cout << "\n";
		}
	}
};*/

/*
struct Node {
	double x, y;	//Koordy nodów ||| x to szerokoœæ, a y to wysokoœæ
	Node(double x1, double y1) {
		x = x1;
		y = y1;
	}

	void printCoords()
	{
		std::cout << "(x: " << x << ",y: " << y << ")";
	}
};

struct Element {
	int id[4];		//Id Nodów w rogach elementu (Nody poiêdzy którymi jest element)
	Element(int id1, int id2, int  id3, int id4) {
		id[0] = id1;
		id[1] = id2;
		id[2] = id3;
		id[3] = id4;
	}

	void printids()
	{
		std::cout << "Id0: " << id[0]+1 << "\nId1: " << id[1]+1 << "\nId2: " << id[2]+1 << "\nId3: " << id[3]+1 << "\n";
	}

	int* returnIdsArrayPointer()
	{
		return id;
	}
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

	Grid(double H1, double B1, int nH1, int nB1) {
		H = H1;
		B = B1;
		nH = nH1;
		nB = nB1;
		nN = nH * nB;
		nE = (nH - 1) * (nB - 1);

		deltaH = H / (nH - 1);
		deltaB = B / (nB - 1);


		for (int i = 0; i < nB; i++) {
			for (int j = 0; j < nH; j++) {
				nodes.push_back(Node(i * deltaB, j * deltaH));
			}
		}
		for (int i = 0; i < (nB - 1); i++) {
			for (int j = 0; j < (nH - 1); j++) {
				int ID1 = (j + (nH)*i);
				int ID2 = ID1 + nH;
				int ID3 = ID2 + 1;
				int ID4 = ID1 + 1;
				elements.push_back(Element(ID1, ID2, ID3, ID4));
			}
		}
	}

	void printAllElements() {
		for (int i = 0; i < elements.size(); i++) {
			std::cout << "Element " << i+1 << ":\n";
			elements[i].printids();
			
		}
	}

	void printEntireGrid() {
		for (int i = 0; i < elements.size(); i++) {
			int* id = elements[i].returnIdsArrayPointer();

			std::cout << "Element " << i + 1 << ":";

			for (int j = 0; j < 4; j++)
			{
				std::cout << "\nId" << j +1 << ": "<< id[j] + 1;
				nodes[id[j]].printCoords();
			}

			std::cout << "\n\n";

		}

	}
};*/

int main() {

	/*Grid test = Grid(0.2,0.1,4,3);
	test.printEntireGrid();
	//std::cout << Gauss1D(1) << "\n";
	//std::cout << Gauss1D(2) << "\n";
	std::cout << Gauss2D(1) << "\n";
	std::cout << Gauss2D(2) << "\n";

	Element4_2D eleme = Element4_2D();
	Element9_2D eleme1 = Element9_2D();
	eleme.printNKsi();
	std::cout << "===========================================================\n";
	eleme.printNEta();
	std::cout << "===========================================================\n";
	std::cout << "======================= Element9_2D =======================\n";
	std::cout << "===========================================================\n";

	eleme1.printNKsi();
	std::cout << "===========================================================\n";
	eleme1.printNEta();*/


	/*std::vector<std::vector<double>> jakob;
	std::vector<std::vector<double>> jakobInversed;*/

	/*double jakob[2][2];
	double jakobInversed[2][2];

	int nPC;
	//Grid grid=Grid(0.025,0.025,2,2);
	Grid grid = Grid(0.05, 0.025, 3, 2);
		nPC = 4;

		for (int i = 0; i < grid.nE; i++) {
			for (int j = 0; j < nPC; j++) {
				if (nPC == 4) {
					jakobian(i, j, jakob, jakobInversed, Element4_2D(), grid);
					
					std::cout << "\n====================================================\n";
					for (int d = 0; d < 2; d++) {
						for (int e = 0; e < 2; e++) {
							std::cout << jakob[d][e] << " \t\t ";
						}
						std::cout << "\n";
					}

				}
				else if (nPC == 9) {
					//jakobian(i, j, jakob, jakobInversed, Element9_2D(), grid);
				}
				else {
					std::cout << "Nope.";
				}

			}
		}*/
	 
	//Sth sth = Sth(Grid(0.025,0.025,2,2));
	Grid(0.025, 0.025, 2, 2);




	return 0;
}