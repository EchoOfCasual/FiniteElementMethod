#include "Elements_2D.h"
#include <math.h>
#include <iostream>
#include "GridPack.h"


Element9_2D::Element9_2D() {
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

void Element9_2D::printNEta() {
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

void Element9_2D::printNKsi() {
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


Element4_2D::Element4_2D() {
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

void Element4_2D::printNEta() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNEta[i][j] << "\t";
			}
			std::cout << "\n";
		}
	}

void Element4_2D::printNKsi() {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				std::cout << derivNKsi[i][j] << "\t";
			}
			std::cout << "\n";
		}
	}