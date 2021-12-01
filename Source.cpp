#include <iostream>
#include "Gauss.h"
#include <vector>
#include "GridPack.h"


int main() {

	//Grid grido = Grid(0.025, 0.025, 2, 2, 4);
	//grido.elements[0].printNEta();
	//grido.elements[0].printHbc();


	//Reminder what stuff in the Grid constructor means (height, width, nodesInHeight, nodesInWidth, integrationPoints|4 or 9|, kFactor|default:30|, aplhaFactor|default:25|, surroundingTemperatureFactor|default:1200|)

	//Grid grido1 = Grid(0.1, 0.1, 4, 4, 4, 25.0, 25.0, 1200.0);	
	Grid grido1 = Grid(0.1, 0.1, 4, 4, 4, 25.0, 300.0, 1200.0, 7800.0, 700.0);
	//grido1.elements[8].printH();
	//grido1.elements[0].printHbc();
	//grido1.elements[0].printNEta();
	//grido1.printHGlobal();

	
	/*grido1.printPGlobal();
	
	for (int i = 0; i < 9; i++) {
		std::cout << " ===== " << i << " element =====";
		grido1.elements[i].printP();

	}*/

	//grido1.solution_t();
	//grido1.printT();
	grido1.printCGlobal();

	return 0;
}