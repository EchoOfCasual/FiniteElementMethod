#pragma once
#include "GridPack.h"

class Element9_2D : public Element {
public:
	double derivNEta[9][4];
	double derivNKsi[9][4];

	Element9_2D();

	void printNEta();
	void printNKsi();
};

class Element4_2D : public Element {
public:
	double derivNEta[4][4];
	double derivNKsi[4][4];

	Element4_2D();

	void printNEta();
	void printNKsi();
};