#if !defined(_FILEPARSER_H_)
#define _FILEPARSER_H_



#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<vector>
#include "vector2.h"
#include "nurbs.h"

using namespace std;

//void readDataFromFile(const char *strFileName, VECTOR2 *&vCP, double *&dKv, int &num);
//
void readNURBSDataFromFile(const char *strFileName, vector<CNurbs*> &curves);
//
//void readPointsFromFile(const char *strFileName, vector<VECTOR2> &testPoints);

void readObjectVertices(const char *sFileName, vector<VECTOR2> &vPoints);

void readVoronoiVertices(const char *sFileName, vector<VECTOR2> &vAllVertices);

void readVoronoiEdges(const char* sFileName);


#endif
