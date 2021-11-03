#include "Gauss.h"
#include <math.h>

double functionForGauss(double x){
	return (5 * x * x + 3 * x + 6);
}

double functionForGaussTwoD(double x, double y) {
	return (5 * x * x * y * y + 3 * x * y + 6);
}

double Gauss1D(int switcher)
{

	if (switcher == 1) {
		double weight[2] = { 1, 1 };
		double nodes[2] = { (-1 / sqrt(3)), (1 / sqrt(3)) };
		return (functionForGauss(nodes[1]) + functionForGauss(nodes[0]));
	}
	else {
		double weight[3] = { 5.0/9.0, 8.0/9.0, 5.0/9.0 };
		double nodes[3] = { -sqrt((3.0/5.0)), 0, sqrt((3.0 / 5.0)) };
		double result = 0;

		for (int i = 0; i < 3; i++) {
			result = result + weight[i] * functionForGauss(nodes[i]);
		}
		
		return result;

	}


};

double Gauss2D(int switcher) {

	if (switcher == 1) {
		double weight[2] = { 1, 1 };
		double nodes[2] = { (-1 / sqrt(3)), (1 / sqrt(3)) };
		double result = 0;

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				result = result + weight[i] * weight[j] * functionForGaussTwoD(nodes[j], nodes[i]);
			}
		}

		return result;
	}
	else {
		double weight[3] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
		double nodes[3] = { -sqrt((3.0 / 5.0)), 0, sqrt((3.0 / 5.0)) };
		double result = 0;

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				result = result + weight[i] * weight[j] * functionForGaussTwoD(nodes[j], nodes[i]);
			}
		}

		return result;

	}

}

