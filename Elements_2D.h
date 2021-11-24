#pragma once
#include "GridPack.h"
#include <iostream>
#include "GridPack.h"

class Element9_2D : public Element {
public:
	double derivNEta[9][4];
	double derivNKsi[9][4];

	Element9_2D(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes);

	void printNEta();
	void printNKsi();
};

class Element4_2D : public Element {
public:
	double derivNEta[4][4];
	double derivNKsi[4][4];

	Element4_2D(int id1, int id2, int  id3, int id4, double k, std::vector<Node> nodes);

	void printNEta();
	void printNKsi();
};