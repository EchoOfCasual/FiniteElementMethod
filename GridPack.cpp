#include "GridPack.h"
#include<vector>
#include<iostream>




Node::Node(double x1, double y1) {
		x = x1;
		y = y1;
	}

void Node::printCoords()
	{
		std::cout << "(x: " << x << ",y: " << y << ")";
	}


	
Element::Element(int id1, int id2, int  id3, int id4) {
		id[0] = id1;
		id[1] = id2;
		id[2] = id3;
		id[3] = id4;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				h[i][j] = 0;
			}
		}
	}

void Element::printids()
	{
		std::cout << "Id0: " << id[0] + 1 << "\nId1: " << id[1] + 1 << "\nId2: " << id[2] + 1 << "\nId3: " << id[3] + 1 << "\n";
	}

int* Element::returnIdsArrayPointer()
	{
		return id;
	}



Grid::Grid(double H1, double B1, int nH1, int nB1, int n, double k) {
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

		//Elementy =====================================================================================================
		double jakob[2][2];
		double jakobInversed[2][2];
		int nPC;


		nPC = n;

		for (int i = 0; i < nE; i++) {
			for (int j = 0; j < nPC; j++) {
				//========

				if (nPC == 4) {
					jakobian(i, j, jakob, jakobInversed, Element4_2D());
				}
				else if (nPC == 9) {
					jakobian(i, j, jakob, jakobInversed, Element9_2D());
				}
				else {
					std::cout << "Nope.";
				}

				double det = jakob[0][0] * jakob[1][1] - jakob[1][0] * jakob[0][1];		//Wyznacznik jakobianu
				double detInv = 1.0 / det;												//Odwrocony wyznacznik jakobianu

				std::cout << "\ndet: " << det << " ||| detInv: " << detInv << "\n============================================\n";
				for (int d = 0; d < 2; d++) {
					for (int e = 0; e < 2; e++) {
						std::cout << jakob[d][e] << " \t\t ";
					}
					std::cout << "\n";
				}
				std::cout << "--------------------------------------------\n";
				double multiplied[2][2];													//Caly ten output jeden

				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						multiplied[i][j] = jakobInversed[i][j] * detInv;
						std::cout << multiplied[i][j] << "\t\t";
					}
					std::cout << "\n";
				}
				std::cout << "\n\n";
				//============================================== H =========================================================
				double derivNX[4][4];
				double derivNY[4][4];

				if (nPC == 4) {
					//	double derivNEta[4][4];
					//  double derivNKsi[4][4];
					Element4_2D elem;
					for (int i1 = 0; i1 < 4; i1++) {
						for (int j1 = 0; j1 < 4; j1++) {
							derivNX[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[0][0] + elem.derivNEta[j][i1] * multiplied[0][1]) * (elem.derivNKsi[j][j1] * multiplied[0][0] + elem.derivNEta[j][j1] * multiplied[0][1]);
							derivNY[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[1][0] + elem.derivNEta[j][i1] * multiplied[1][1]) * (elem.derivNKsi[j][j1] * multiplied[1][0] + elem.derivNEta[j][j1] * multiplied[1][1]);
						}
					}
				}
				else if (nPC == 9) {
					Element9_2D elem;
					for (int i1 = 0; i1 < 4; i1++) {
						for (int j1 = 0; j1 < 4; j1++) {
							derivNX[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[0][0] + elem.derivNEta[j][i1] * multiplied[0][1]) * (elem.derivNKsi[j][j1] * multiplied[0][0] + elem.derivNEta[j][j1] * multiplied[0][1]);
							derivNY[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[1][0] + elem.derivNEta[j][i1] * multiplied[1][1]) * (elem.derivNKsi[j][j1] * multiplied[1][0] + elem.derivNEta[j][j1] * multiplied[1][1]);
						}
					}
				}
				else {
					std::cout << "Nope.";
				}

				double hpc[4][4];

				for (int i1 = 0; i1 < 4; i1++) {
					for (int j1 = 0; j1 < 4; j1++) {
						hpc[i1][j1] = (derivNX[i1][j1] + derivNY[i1][j1]) * k * det;
						std::cout << hpc[i1][j1] << "\t";
						elements[i].h[i1][j1] += hpc[i1][j1];
					}
					std::cout << "\n";
				}





			}
			std::cout << "\n====== H =======\n";
			for (int i1 = 0; i1 < 4; i1++) {
				for (int j1 = 0; j1 < 4; j1++) {
					std::cout << elements[i].h[i1][j1] << "\t";
				}
				std::cout << "\n";
			}
			std::cout << "\n====================================== Next Element ==================================\n";
		}

	}

void Grid::printAllElements() {
		for (int i = 0; i < elements.size(); i++) {
			std::cout << "Element " << i + 1 << ":\n";
			elements[i].printids();

		}
	}

void Grid::printEntireGrid() {
		for (int i = 0; i < elements.size(); i++) {
			int* id = elements[i].returnIdsArrayPointer();

			std::cout << "Element " << i + 1 << ":";

			for (int j = 0; j < 4; j++)
			{
				std::cout << "\nId" << j + 1 << ": " << id[j] + 1;
				nodes[id[j]].printCoords();
			}

			std::cout << "\n\n";

		}

	}

void Grid::jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element9_2D elem) {
	//The very same imp belllowww, not tested yet tho, soooooooo, good luck with that future me
	jakobian[0][0] = 0;
	jakobian[0][1] = 0;
	jakobian[1][0] = 0;
	jakobian[1][1] = 0;

	//double xHardcode[4] = { 0, 0.025, 0.025, 0 };
	//double yHardcode[4] = { 0, 0, 0.025, 0.025 };

	for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].y;
		jakobian[1][0] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobian[1][1] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].y;
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
		jakobianInv[1][1] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobianInv[0][1] -= elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].y;
		jakobianInv[1][0] -= elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobianInv[0][0] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].y;
	}

};

void Grid::jakobian(int elemI, int pCJ, /*std::vector<std::vector<double>> &jakobian*/ double jakobian[][2], /* std::vector<std::vector<double>>& jakobianInv*/ double jakobianInv[][2], Element4_2D elem) {
	jakobian[0][0] = 0;
	jakobian[0][1] = 0;
	jakobian[1][0] = 0;
	jakobian[1][1] = 0;

	//double xHardcode[4] = { 0, 0.025, 0.025, 0 };
	//double yHardcode[4] = { 0, 0, 0.025, 0.025 };

	for (int i = 0; i < 4; i++) {
		jakobian[0][0] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobian[0][1] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].y;
		jakobian[1][0] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobian[1][1] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].y;
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
		jakobianInv[1][1] += elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobianInv[0][1] -= elem.derivNKsi[pCJ][i] * nodes[elements[elemI].id[i]].y;
		jakobianInv[1][0] -= elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].x;
		jakobianInv[0][0] += elem.derivNEta[pCJ][i] * nodes[elements[elemI].id[i]].y;
	}


};

void Grid::calculate(int n, double k) {

	double jakob[2][2];
	double jakobInversed[2][2];
	int nPC;


		nPC = n;

		for (int i = 0; i < nE; i++) {
			for (int j = 0; j < nPC; j++) {
				//========

				if (nPC == 4) {
					jakobian(i, j, jakob, jakobInversed, Element4_2D());
				}
				else if (nPC == 9) {
					jakobian(i, j, jakob, jakobInversed, Element9_2D());
				}
				else {
					std::cout << "Nope.";
				}

				double det = jakob[0][0] * jakob[1][1] - jakob[1][0] * jakob[0][1];		//Wyznacznik jakobianu
				double detInv = 1.0 / det;												//Odwrocony wyznacznik jakobianu

				std::cout << "\ndet: " << det << " ||| detInv: " << detInv << "\n============================================\n";
				for (int d = 0; d < 2; d++) {
					for (int e = 0; e < 2; e++) {
						std::cout << jakob[d][e] << " \t\t ";
					}
					std::cout << "\n";
				}
				std::cout << "--------------------------------------------\n";
				double multiplied[2][2];													//Caly ten output jeden

				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						multiplied[i][j] = jakobInversed[i][j] * detInv;
						std::cout << multiplied[i][j] << "\t\t";
					}
					std::cout << "\n";
				}
				std::cout << "\n\n";
				//============================================== H =========================================================
				double derivNX[4][4];
				double derivNY[4][4];

				if (nPC == 4) {
					//	double derivNEta[4][4];
					//  double derivNKsi[4][4];
					Element4_2D elem;
					for (int i1=0; i1 < 4; i1++) {
						for (int j1=0; j1 < 4; j1++) {
							derivNX[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[0][0] + elem.derivNEta[j][i1] * multiplied[0][1]) * (elem.derivNKsi[j][j1] * multiplied[0][0] + elem.derivNEta[j][j1] * multiplied[0][1]);
							derivNY[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[1][0] + elem.derivNEta[j][i1] * multiplied[1][1]) * (elem.derivNKsi[j][j1] * multiplied[1][0] + elem.derivNEta[j][j1] * multiplied[1][1]);
						}
					}
				}
				else if (nPC == 9) {
					Element9_2D elem;
					for (int i1=0; i1 < 4; i1++) {
						for (int j1=0; j1 < 4; j1++) {
							derivNX[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[0][0] + elem.derivNEta[j][i1] * multiplied[0][1]) * (elem.derivNKsi[j][j1] * multiplied[0][0] + elem.derivNEta[j][j1] * multiplied[0][1]);
							derivNY[i1][j1] = (elem.derivNKsi[j][i1] * multiplied[1][0] + elem.derivNEta[j][i1] * multiplied[1][1]) * (elem.derivNKsi[j][j1] * multiplied[1][0] + elem.derivNEta[j][j1] * multiplied[1][1]);
						}
					}
				}
				else {
					std::cout << "Nope.";
				}

				double hpc[4][4];

				for (int i1=0; i1 < 4; i1++) {
					for (int j1=0; j1 < 4; j1++) {
						hpc[i1][j1] = (derivNX[i1][j1] + derivNY[i1][j1])* k  *det;
						std::cout << hpc[i1][j1] << "\t";
						elements[i].h[i1][j1] += hpc[i1][j1];
					}
					std::cout << "\n";
				}





			}
			std::cout << "\n====== H =======\n";
			for (int i1 = 0; i1 < 4; i1++) {
				for (int j1 = 0; j1 < 4; j1++) {
					std::cout << elements[i].h[i1][j1] << "\t";
				}
				std::cout << "\n";
			}
			std::cout << "\n====================================== Next Element ==================================\n";
		}

}