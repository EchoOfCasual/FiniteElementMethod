#include <iostream>
#include "Gauss.h"
#include <vector>
#include "GridPack.h"


int main() {

	//Grid grido = Grid(0.025, 0.025, 2, 2, 4);
	//grido.elements[0].printNEta();
	//grido.elements[0].printHbc();



	//Grid grido1 = Grid(0.1, 0.1, 4, 4, 4, 25.0, 25.0, 1200.0);
	Grid grido1 = Grid(0.1, 0.1, 4, 4, 9, 25.0, 25.0, 1200.0);
	//grido1.elements[0].printHbc();
	//grido1.elements[0].printNEta();
	//grido1.printHGlobal();

	grido1.printPGlobal();
	
	for (int i = 0; i < 9; i++) {
		std::cout << " ===== " << i << " element =====";
		grido1.elements[i].printP();

	}

	return 0;
}