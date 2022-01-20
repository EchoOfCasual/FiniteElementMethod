#include "GridPack.h"
#include<vector>
#include<iostream>
#include<math.h>

#include <algorithm>



Node::Node(double x1, double y1, short int bc1) {
		x = x1;						//Coord x
		y = y1;						//Coord y
		bc = bc1;					//Does it have a border condition flag
	}

void Node::print()
	{
		std::cout << "(x: " << x << ",y: " << y << ") Border condition: "<<bc;
	}


Element::Element(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, int points, double(*derivNEta1)[4], double(*derivNKsi1)[4], double alpha1, double hbc1[4][4], double tSurrounding1, double p1[4], double ro1, double specHeat1, double(*N1)[4]) {
	id[0] = id1;					//Initialising ids of nodes creating element
	id[1] = id2;
	id[2] = id3;
	id[3] = id4;
	nIntegrationPoints = points;	//How many integration points there is (4 or 9 in 2d)



	derivNEta = new double[nIntegrationPoints][4];
	derivNKsi = new double[nIntegrationPoints][4];
	N = new double[nIntegrationPoints][4];

	//derivNEta = derivNEta1;			//Assigning arrays of deriv N/Eta and N/Ksi provied by static funtion
	//derivNKsi = derivNKsi1;			//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Ogarnij, bo sie nie wpisuje do elementu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Juz zapisuje


	for (int i1 = 0; i1 < nIntegrationPoints; i1++) { //Assigning arrays of deriv N/Eta and N/Ksi and N (shape functions) provied by static funtion, line 31 does not work
		for (int j1 = 0; j1 < 4; j1++) {
			derivNEta[i1][j1] = derivNEta1[i1][j1];
			derivNKsi[i1][j1] = derivNKsi1[i1][j1];
			N[i1][j1] = N1[i1][j1];
		}
	}

	k1 = k;											// Assigning valuses k and alpha from static...
	alpha = alpha1;

	for (int i1 = 0; i1 < 4; i1++) {				// Assigning valuses p and hbc from static & initialising h & c
		p[i1] = p1[i1];
		for (int j1 = 0; j1 < 4; j1++) {
			h[i1][j1] = 0;
			c[i1][j1] = 0;
			hbc[i1][j1] = hbc1[i1][j1];
		}
	}

	tSurrounding = tSurrounding1; //Setting temperature of the sourrounding area
	ro = ro1;					// Assigning ro form static and specific heat
	specHeat = specHeat1;

	// ===================================== Calculation of H ======================================
	for (int j = 0; j < nIntegrationPoints; j++) {
		//======== Calculating jacobian (and jacobian inversed) =============
		jacobian[0][0] = 0;			//Setting jakobian matrix to 0
		jacobian[0][1] = 0;
		jacobian[1][0] = 0;
		jacobian[1][1] = 0;
		
									//Calculating jacobian (x and y for both xsi and eta) in 4 nodes
		for (int i = 0; i < 4; i++) {
			jacobian[0][0] += derivNKsi[j][i] * nodes[id[i]].x;
			jacobian[0][1] += derivNKsi[j][i] * nodes[id[i]].y;
			jacobian[1][0] += derivNEta[j][i] * nodes[id[i]].x;
			jacobian[1][1] += derivNEta[j][i] * nodes[id[i]].y;
		}


									//Setting jakobian inversed matrix to 0
		jacobianInversed[0][0] = 0;
		jacobianInversed[0][1] = 0;
		jacobianInversed[1][0] = 0;
		jacobianInversed[1][1] = 0;
				

									//Calculating inversed jacobian (note to myself: consider using calculated jacobian earlier to calculate this - it is more optimized tho i neeed to sleep)
		for (int i = 0; i < 4; i++) {
			jacobianInversed[1][1] += derivNKsi[j][i] * nodes[id[i]].x;
			jacobianInversed[0][1] -= derivNKsi[j][i] * nodes[id[i]].y;
			jacobianInversed[1][0] -= derivNEta[j][i] * nodes[id[i]].x;
			jacobianInversed[0][0] += derivNEta[j][i] * nodes[id[i]].y;
		}
		//=====================================================================


		jacobDet = jacobian[0][0] * jacobian[1][1] - jacobian[1][0] * jacobian[0][1];	//Det of jacobian
		jacobDetInversed = 1.0 / jacobDet;												//Inversed det of jacobian
																						//Just printing
		/*std::cout << "\ndet: " << jacobDet << " ||| detInv: " << jacobDetInversed << "\n============================================\n";
		for (int d = 0; d < 2; d++) {
			for (int e = 0; e < 2; e++) {
				std::cout << jacobian[d][e] << " \t\t ";
			}
			std::cout << "\n";
		}
		std::cout << "--------------------------------------------\n";*/
		double multiplied[2][2];													//Inversed jacobian multiplied by inversed det of jacobian

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				multiplied[i][j] = jacobianInversed[i][j] * jacobDetInversed;		//Above
				//std::cout << multiplied[i][j] << "\t\t";							//Aaaaand some printing as well
			}
			//std::cout << "\n";
		}
		//std::cout << "\n\n";

		//============================================== Calculating H (but first deriv N/x and N/y) =========================================================
		double derivNX[4][4];		//Helpful matrixes (1 integration point, vector multiplied by trnapsosed itself)
		double derivNY[4][4];
		double Ntemp[4][4];

		for (int i1 = 0; i1 < 4; i1++) {
			for (int j1 = 0; j1 < 4; j1++) {	//It is actually derivNX * derivNXW^T and the same with Y, just a note
				derivNX[i1][j1] = (derivNKsi[j][i1] * multiplied[0][0] + derivNEta[j][i1] * multiplied[0][1]) * (derivNKsi[j][j1] * multiplied[0][0] + derivNEta[j][j1] * multiplied[0][1]);	//Exactly what already written above. Note: multiplied[0] is for x and multiplied[1] is for y. More visable on PDF
				derivNY[i1][j1] = (derivNKsi[j][i1] * multiplied[1][0] + derivNEta[j][i1] * multiplied[1][1]) * (derivNKsi[j][j1] * multiplied[1][0] + derivNEta[j][j1] * multiplied[1][1]);

				Ntemp[i1][j1] = N[j][i1] * N[j][j1]; //Needed to get c
				//std::cout << Ntemp[i1][j1] << "\t";
			}
				//std::cout << "\n";
		}
		//std::cout << "\n============================\n";

		double weight = 1;
		if (nIntegrationPoints == 9) {
			if (j == 0 || j == 2 || j == 6 || j == 8) {
				weight = 25.0 / 81.0;
			}
			else if(j == 4) {
				weight = 64.0 / 81.0;
			}
			else
			{
				weight = 40.0 / 81.0;
			}
		}

		double hpc[4][4];			//Well, h but in 1 integration point!


		//std::cout << "\n\n==============\n";
		for (int i1 = 0; i1 < 4; i1++) {
			for (int j1 = 0; j1 < 4; j1++) {
				hpc[i1][j1] = weight * (derivNX[i1][j1] + derivNY[i1][j1]) * k * jacobDet;	//Comment above.
				//std::cout << hpc[i1][j1] << "\t";	//Some printing here as well
				h[i1][j1] += hpc[i1][j1];	//One H to rule them all, One H to find them, One H to bring them all and in the element bind them. (In one Element - isnt as powerful as hGlobal)


				//std::cout << "weight: " << weight <<" || ro: " << ro <<" || specHeat: " << specHeat <<" || Ntemp[i1][j1]: " << Ntemp[i1][j1] <<"\n";
				c[i1][j1] += weight * ro * specHeat * Ntemp[i1][j1] * jacobDet;
			}
			//std::cout << "\n";
		}


	}
	/*std::cout << "\n====== H =======\n";	//Printing the one H
	for (int i1 = 0; i1 < 4; i1++) {
		for (int j1 = 0; j1 < 4; j1++) {
			std::cout << h[i1][j1] << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n====================================== Next Element ==================================\n";*/
}

/*Element::~Element() {
	delete derivNEta;
	delete derivNKsi;
}*/

Element* Element::create9(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, double alpha1, double tSurrounding1, double ro1, double specHeat1) {		//This thing in here is to create element9_2D, some prep for derivs N/eta and N/xsi and then invoking the actuall constructor

	double derivNEta1[9][4];
	double derivNKsi1[9][4];
	double N1[9][4];

	for (int i = 0; i < 9; i++) {

		double eta, ksi;

		if (i == 0 || i == 3 || i == 6) {
			ksi = -sqrt(3.0 / 5.0);
		}
		else if (i == 1 || i == 4 || i == 7) {
			ksi = 0;
		}
		else {
			ksi = sqrt(3.0 / 5.0);
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

		derivNEta1[i][0] = -0.25 * (1 - ksi);
		derivNEta1[i][1] = -0.25 * (1 + ksi);
		derivNEta1[i][2] = 0.25 * (1 + ksi);
		derivNEta1[i][3] = 0.25 * (1 - ksi);

		derivNKsi1[i][0] = -0.25 * (1 - eta);
		derivNKsi1[i][1] = 0.25 * (1 - eta);
		derivNKsi1[i][2] = 0.25 * (1 + eta);
		derivNKsi1[i][3] = -0.25 * (1 + eta);

		N1[i][0] = 0.25 * (1 - eta) * (1 - ksi);
		N1[i][1] = 0.25 * (1 - eta) * (1 + ksi);
		N1[i][2] = 0.25 * (1 + eta) * (1 + ksi);
		N1[i][3] = 0.25 * (1 + eta) * (1 - ksi);

	}

	// ================================== Calculation of hbc and p ========================================= Czy przeniesc do wspolnej funkcji??? I teraz tez P (nope: niby kilka powtorzonych linijek, ale nie ma czasu i bedzie sie gorzej czytac)
	double hbc1[4][4] = { {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0} };
	double p1[4] = {0,0,0,0};
	int ids[4] = { id1, id2, id3, id4 };

	for (int i = 0; i < 4; i++) {
		if (nodes[ids[i]].bc == 1 && nodes[ids[(i + 1) % 4]].bc == 1) {
			double tempPc1[4] = { 0,0,0,0 };
			double tempPc2[4] = { 0,0,0,0 };
			double tempPc3[4] = { 0,0,0,0 };

			tempPc1[i] = (1.0 - sqrt(3.0/5.0)) / 2.0;			//Cakculating x of point in 0 -> 1 (instead of -1 -> 1) (from -1 -> 1 Gauss to 0 -> 1 ez representation of shape func value)
			tempPc1[(i + 1) % 4] = 1 - tempPc1[i];

			tempPc2[(i + 1) % 4] = 0.5;
			tempPc2[i] = 0.5;

			tempPc3[(i + 1) % 4] = tempPc1[i];
			tempPc3[i] = tempPc1[(i + 1) % 4];

			double detJ1 = (sqrt(pow(nodes[ids[i]].x - nodes[ids[(i + 1) % 4]].x, 2) + pow(nodes[ids[i]].y - nodes[ids[(i + 1) % 4]].y, 2))) / 2.0;	//Dlugosc boku przez 2 det jakobianu || (length of the element side divided by 2) = (Jacobian det)

			for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
				double weight = 5.0 / 9.0;
				p1[hbcRow] += weight * tempPc1[hbcRow]* tSurrounding1 * alpha1 * detJ1;
				p1[hbcRow] += (8.0 / 9.0) * tempPc2[hbcRow] * tSurrounding1 * alpha1 * detJ1;
				p1[hbcRow] += weight * tempPc3[hbcRow] * tSurrounding1 * alpha1 * detJ1;
				for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
					
					hbc1[hbcRow][hbcColumn] += weight * tempPc1[hbcRow] * tempPc1[hbcColumn] * alpha1 * detJ1;
					hbc1[hbcRow][hbcColumn] += (8.0 / 9.0) * tempPc2[hbcRow] * tempPc2[hbcColumn] * alpha1 * detJ1; // Weight in middle integration point * matrix * transposed matrix *alpha * det
					hbc1[hbcRow][hbcColumn] += weight * tempPc3[hbcRow] * tempPc3[hbcColumn] * alpha1 * detJ1;	// Weight in the outer integration point * matrix * transposed matrix *alpha * det
				}
			}
		}
	}
	/*for (int i1 = 0; i1 < 4; i1++) {
		for (int j1 = 0; j1 < 4; j1++) {
			std::cout << hbc1[i1][j1] << "\t";
		}
		std::cout << "\n";
	}*/
	// ===============================================================================================
	return new Element(id1, id2, id3, id4, k, nodes, 9, derivNEta1, derivNKsi1, alpha1, hbc1, tSurrounding1, p1, ro1, specHeat1, N1);	//Invoking and returning constructor
}

Element* Element::create4(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes, double alpha1, double tSurrounding1, double ro1, double specHeat1) {		//This thing in here is to create element4_2D, some prep for derivs N/eta and N/xsi and then invoking the actuall constructor

	double derivNEta1[4][4];
	double derivNKsi1[4][4];
	double N1[4][4];


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

		derivNEta1[i][0] = -0.25 * (1 - ksi);		// deriv N/eta in integr points (0.25 not 0.5 cause in 0 -> 1)
		derivNEta1[i][1] = -0.25 * (1 + ksi);
		derivNEta1[i][2] = 0.25 * (1 + ksi);
		derivNEta1[i][3] = 0.25 * (1 - ksi);

		derivNKsi1[i][0] = -0.25 * (1 - eta);
		derivNKsi1[i][1] = 0.25 * (1 - eta);
		derivNKsi1[i][2] = 0.25 * (1 + eta);
		derivNKsi1[i][3] = -0.25 * (1 + eta);


		N1[i][0] = 0.25 * (1 - eta) * (1 - ksi);	//Shape func in integration points
		N1[i][1] = 0.25 * (1 - eta) * (1 + ksi);
		N1[i][2] = 0.25 * (1 + eta) * (1 + ksi);
		N1[i][3] = 0.25 * (1 + eta) * (1 - ksi);
	}

	// ================================== Calculation of hbc & p =========================================
	double hbc1[4][4] = { {0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0} };
	double p1[4] = { 0,0,0,0 };
	int ids[4] = { id1, id2, id3, id4 };

	for (int i = 0; i < 4; i++) {
		if (nodes[ids[i]].bc == 1 && nodes[ids[(i + 1) % 4]].bc == 1) {
			double tempPc1[4] = { 0,0,0,0 };
			double tempPc2[4] = { 0,0,0,0 };

			tempPc1[i] = (1.0 - 1.0 / sqrt(3)) / 2.0;				
			tempPc1[(i + 1) % 4] = 1 - ((1.0 - 1.0 / sqrt(3)) / 2.0);

			tempPc2[(i + 1) % 4] = (1.0 - 1.0 / sqrt(3)) / 2.0;
			tempPc2[i] = 1 - ((1.0 - 1.0 / sqrt(3)) / 2.0);

			double detJ1 = (sqrt( pow(nodes[ids[i]].x - nodes[ids[(i + 1) % 4]].x, 2) + pow(nodes[ids[i]].y - nodes[ids[(i + 1) % 4]].y, 2) )) / 2.0;	//Dlugosc boku przez 2 -> det jakobianu

			for (int hbcRow = 0; hbcRow < 4; hbcRow++) {
				p1[hbcRow] += tempPc1[hbcRow] * tSurrounding1 * alpha1 * detJ1;
				p1[hbcRow] += tempPc2[hbcRow] * tSurrounding1 * alpha1 * detJ1;
				for (int hbcColumn = 0; hbcColumn < 4; hbcColumn++) {
					hbc1[hbcRow][hbcColumn] += tempPc1[hbcRow] * tempPc1[hbcColumn] * alpha1 * detJ1;
					hbc1[hbcRow][hbcColumn] += tempPc2[hbcRow] * tempPc2[hbcColumn] * alpha1 * detJ1;
				}
			}
		}
	}



	return new Element(id1, id2, id3, id4, k, nodes, 4, derivNEta1, derivNKsi1, alpha1, hbc1, tSurrounding1, p1, ro1, specHeat1, N1);	//Invoking and returning constructor
}


void Element::printids()
	{
		std::cout << "Id0: " << id[0] + 1 << "\nId1: " << id[1] + 1 << "\nId2: " << id[2] + 1 << "\nId3: " << id[3] + 1 << "\n";
	}

int* Element::returnIdsArrayPointer()
	{
		return id;
	}


void Element::printNEta() {
	for (int i = 0; i < nIntegrationPoints; i++) {
		//double* temp = derivNEta[i];
		for (int j = 0; j < 4; j++) {
			std::cout << derivNEta[i][j] << "\t";
			//std::cout << temp[j] << "\t";
			//if (i == 1 || i == 4 || i == 7) {
			//	std::cout << "\t";
			//}
		}
		std::cout << "\n";
	}
}

void Element::printNKsi() {
	for (int i = 0; i < nIntegrationPoints; i++) {
		for (int j = 0; j < 4; j++) {
			std::cout << derivNKsi[i][j] << "\t";

			//if (i > 2 && i < 6) {
			//	std::cout << "\t";
			//}

		}
		std::cout << "\n";
	}
}
void Element::printH()
{
	std::cout <<"\n ========================================= H ========================================= \n";
	for (int i1 = 0; i1 < 4; i1++) {
		for (int j1 = 0; j1 < 4; j1++) {
			std::cout << h[i1][j1] << "\t";
		}
		std::cout << "\n";
	}
}
void Element::printHbc()
{
	std::cout << "\n ========================================= Hbc ========================================= \n";
	for (int i1 = 0; i1 < 4; i1++) {
		for (int j1 = 0; j1 < 4; j1++) {
			std::cout << hbc[i1][j1] << "\t";
		}
		std::cout <<"\n";
	}
}

void Element::printP()		// ========================================= P =========================================
{
	std::cout << "\n P: \n";
	for (int i1 = 0; i1 < 4; i1++) {
		std::cout << " " << p[i1] << "\n";
	}
}

void Element::printC()		// ========================================= C =========================================
{
	std::cout << "\n ========================================= C ========================================= \n";
	for (int i1 = 0; i1 < 4; i1++) {
		for (int j1 = 0; j1 < 4; j1++) {
			std::cout << c[i1][j1] << "\t";
		}
		std::cout << "\n";
	}
}




Grid::Grid(double H1, double B1, int nH1, int nB1, int n, double k, double alpha1, double tSurrounding1, double ro1, double specHeat1) {
		H = H1;						//Height
		B = B1;						//Width
		nH = nH1;					//Amount of nodes in height
		nB = nB1;					//Amount of nodes in width
		nN = nH * nB;				//Amount of nodes overall
		nE = (nH - 1) * (nB - 1);	//Amount of elements
		nIntegrationPoints = n;
		k1 = k;
		alpha = alpha1;
		tSurrounding = tSurrounding1;

		deltaH = H / (nH - 1);		//Height of one element (if elements are the same)
		deltaB = B / (nB - 1);		//Width of one element


		for (int i = 0; i < nB; i++) {			//Creating nodes (For the same elements (equal in size)). I mean this entire constructor is based on that assumption soo
			for (int j = 0; j < nH; j++) {
				short int bc = 0;
				if (i == 0 || i == nB-1 || j == 0 || j == nH-1)
				{
					bc = 1;
				}
				nodes.push_back(Node(i * deltaB, j * deltaH, bc));
			}
		}
		for (int i = 0; i < (nB - 1); i++) {	//Creating elements 
			for (int j = 0; j < (nH - 1); j++) {
				int ID1 = (j + (nH)*i);			//Calculating which nodes element is consisted of
				int ID2 = ID1 + nH;				
				int ID3 = ID2 + 1;
				int ID4 = ID1 + 1;
				if (n == 4) {					//Normal case (4 4 elements in 2 d)
					elements.push_back(*Element::create4(ID1, ID2, ID3, ID4, k, /*(not so )temp solution(it seems)*/nodes, alpha1, tSurrounding1, ro1, specHeat1));
				}
				else if (n == 9) {				//Abnormal case (4 9 elements in 2 d) Not working for some devilish reason (it works now)
					elements.push_back(*Element::create9(ID1, ID2, ID3, ID4, k, /*(not so )temp solution(it seems)*/nodes, alpha1, tSurrounding1, ro1, specHeat1));
				}
				else {							//Yeee, nope
					std::cout << "Wrong n! (Integration points)";
					exit(1);
				}
			}
		}

		//============================== calculating h global & p global ============================
		for (int i = 0; i < nN; i++) {							//Initiating hGlobal with 0
			std::vector<double> tempForTheColumns;
			std::vector<double> tempForTheColumnsC;
			pGlobal.push_back(0.0);								//Initialising pGlobal with 0s
			for (int j = 0; j < nN; j++) {
				tempForTheColumns.push_back(0.0);
				tempForTheColumnsC.push_back(0.0);
			}
			hGlobal.push_back(tempForTheColumns);
			cGlobal.push_back(tempForTheColumnsC);
		}

		for (int elementsColumns = 0; elementsColumns < nB - 1; elementsColumns++) {							//setting proper values for hGlobal

			for (int elementRows = 0; elementRows < nH - 1; elementRows++) {
				
				for (int hColumns = 0; hColumns < 4; hColumns++) {							

					for (int hRows = 0; hRows < 4; hRows++) {

						hGlobal[elements[elementsColumns * (nH - 1) + elementRows].id[hColumns]][elements[elementsColumns * (nH - 1) + elementRows].id[hRows]] += elements[elementsColumns * (nH - 1) + elementRows].h[hColumns][hRows] + elements[elementsColumns * (nH - 1) + elementRows].hbc[hColumns][hRows];

						cGlobal[elements[elementsColumns * (nH - 1) + elementRows].id[hColumns]][elements[elementsColumns * (nH - 1) + elementRows].id[hRows]] += elements[elementsColumns * (nH - 1) + elementRows].c[hColumns][hRows];


					}

				}


			}
			
		}

		for (int elementsI = 0; elementsI < nE; elementsI++) {													// setting proper values for pGlobal (Probably better solution than the above)	
			for (int pRows = 0; pRows < 4; pRows++) {
				pGlobal[elements[elementsI].id[pRows]] += elements[elementsI].p[pRows];
			}
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
				nodes[id[j]].print();
			}

			std::cout << "\n\n";

		}

	}

void Grid::printHGlobal() {
	for (int i = 0; i < nN; i++) {

		for (int j = 0; j < nN; j++)
		{
			std::cout << hGlobal[i][j] << " || ";
		}

		std::cout << "\n";

	}

}

void Grid::printPGlobal() {
	std::cout << " ============================ pGlobal ============================ \n";

	for (int i = 0; i < nN; i++) {

		std::cout << " " << i+1 << ": "<< pGlobal[i] << "\n";

	}

}


std::vector<double> Grid::solution_t(double tau, std::vector<double> initialTemp) {
	int i, j, k;
	double m, s;
	int n = hGlobal.size();

	std::vector<double> t;
	t.resize(n, 0);

	std::vector<std::vector<double>> array = hGlobal;

	std::vector<double> initialTempVector;
	initialTempVector.resize(n,0);

	for(int i =0; i < hGlobal.size(); i++) {
		for(int j =0; j < hGlobal.size(); j++) {
			array[i][j] += cGlobal[i][j]/tau;
			initialTempVector[i] += cGlobal[i][j]/tau *initialTemp[j];
		}
		//std::cout << pGlobal[i] + initialTempVector[i] << std::endl;
	}

	//array.push_back(pGlobal);
	for (int i = 0; i < n; i++) {
		array[i].push_back((pGlobal[i]+initialTempVector[i]));
	}

	// eliminacja wsp�czynnik�w
	for (i = 0; i < n - 1; i++)
	{
		for (j = i + 1; j < n; j++)
		{
			if (array[i][i] == 0.0) { std::cout << "Operacja powiod�a sie"; return {}; }
			m = -array[j][i] / array[i][i];
			for (k = i + 1; k <= n; k++)
				array[j][k] += m * array[i][k];
		}
	}

	// wyliczanie niewiadomych

	for (i = n - 1; i >= 0; i--)
	{
		s = array[i][n];
		for (j = n - 1; j >= i + 1; j--)
			s -= array[i][j] * t[j];
		if (array[i][i] == 0.0) { std::cout << "Operacja nie powiod�a sie"; return {}; }
		t[i] = s / array[i][i];
	}

	//std::cout << "Operacja powiod�a sie";
	return t;

}

std::vector<double> Grid::final_solution_t(double tau, std::vector<double> initialTemp, int iterationsNumber) {


	std::vector<double> output = initialTemp;

	for (int i = 0; i < iterationsNumber; i++) {
		/*Grid newValuesCalculation = Grid(H, B, nH, nB, nIntegrationPoints, k1, alpha, tSurrounding, 7800.0, 700.0); // w Tym przypadku pozostaj� takie same, ale og�lnie g�sto�� i ciep�o w�a�ciwe mog� si� zmienia�. W tym miejscu okazuje sie, �e lepszym wzorcem od metody wytw�rczej by�by builder, wcze�niej nie wiedzia�em, �e b�dzie wykorzystywane co� takiego, a teraz niestety refaktoryzacja zaj�aby za du�o czasu
		nodes = newValuesCalculation.nodes;
		elements = newValuesCalculation.elements;
		hGlobal = newValuesCalculation.hGlobal;
		cGlobal = newValuesCalculation.cGlobal;
		pGlobal = newValuesCalculation.pGlobal;*/

		output = solution_t(tau, output);

		std::cout << "======== " << i + 1 << " ========" << std::endl;
		std::cout << "\nMin Element = "
			<< *std::min_element(output.begin(), output.end());


		std::cout << "\nMax Element = "
			<< *std::max_element(output.begin(), output.end());
		
		/*std::cout << "======== " << i+1 << " ========" << std::endl;
		for (int i = 0; i < hGlobal.size(); i++) {
			std::cout << output[i] << std::endl;
		}*/
		

	}
	return output;

}

/*void Grid::printT() {
	std::cout << " ============================ t ============================ \n";

	for (int i = 0; i < nN; i++) {

		std::cout << " " << i + 1 << ": " << t[i] << "\n";

	}

}*/

void Grid::printCGlobal() {
	for (int i = 0; i < nN; i++) {

		for (int j = 0; j < nN; j++)
		{
			std::cout << cGlobal[i][j] << " || ";
		}

		std::cout << "\n";

	}

}