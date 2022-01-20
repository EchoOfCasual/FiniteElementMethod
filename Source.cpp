#include <iostream>
#include "Gauss.h"
#include <vector>
#include "GridPack.h"


int main() {

	//Grid grido = Grid(0.025, 0.025, 2, 2, 4);
	//grido.elements[0].printNEta();
	//grido.elements[0].printHbc();


	//Reminder what stuff in the Grid constructor means (height, width, nodesInHeight, nodesInWidth, integrationPoints|4 or 9|, kFactor|default:30|, aplhaFactor|default:25|, surroundingTemperatureFactor|default:1200|, ro(gestosc), specHeat(cieplo wlasciwe))

	//Grid grido1 = Grid(0.1, 0.1, 4, 4, 4, 25.0, 25.0, 1200.0);	
	Grid grido1 = Grid(0.1, 0.1, 31, 31, 9, 25.0, 300.0, 1200.0, 7800.0, 700.0);
	//grido1.elements[8].printH();
	//grido1.elements[0].printHbc();
	//grido1.elements[0].printNEta();
	//grido1.printHGlobal();
	//std::cout << "========================================================================================\n";
	//grido1.printCGlobal();
	std::vector<double> init_temp;
	init_temp.resize(grido1.hGlobal.size(), 100.0);

	std::vector<double> output = grido1.final_solution_t(1, init_temp, 20);

	//grido1.printEntireGrid();
	/*grido1.printPGlobal();
	
	for (int i = 0; i < 9; i++) {
		std::cout << " ===== " << i << " element =====";
		grido1.elements[i].printP();

	}*/

	//grido1.solution_t();
	//grido1.printT();
	//grido1.printHGlobal();
	//grido1.elements[3].printC();


	//for (int i = 0; i < grido1.hGlobal.size(); i++) {
	//	std::cout << output[i] << std::endl;
	//}

	//grido1.printT();
	

	return 0;
}