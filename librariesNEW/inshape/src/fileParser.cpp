//
//
// functions to read BSpline data from a file
//
//
//
//
#include <fileParser.h>


//void readDataFromFile(const char *strFileName, VECTOR2 *&vCp, double *&dKv, int &num)
//{
//
//	int i;
//
//	/* open stream for reading */
//	ifstream file(strFileName);
//
//	if(!(file.is_open()))
//	{
//		cerr<<"Could not open file "<<strFileName<<endl;
//		return;
//	}
//
//	char temp[100];
//	char *pEnd;
//	double dbl;
//
//	/* read number of points */
//	file.getline(temp, 100);
//	int numPoints = atoi(temp);
//
//	/* assign num value in main program */
//	num = numPoints;
//
//	/* allocate memory for array of control points */
//	vCp = new VECTOR2[numPoints]; 
//	cout<<"reading points"<<endl;
//	/* read in control points one by one */
//	for(i = 0; i < numPoints; i++)
//	{
//		file.getline(temp, 100);
//		dbl=strtod(temp, &pEnd);
//		vCp[i].x = dbl;
//		file.getline(temp, 100);
//		dbl=strtod(temp, &pEnd);
//		vCp[i].y = dbl;
//		cout<<vCp[i];
//		
//	}//end for i
//
//	/* read number of knots */
//	file.getline(temp, 100);
//	int numKnots = atoi(temp);
//
//	/* allocate memory for knot vector */
//	dKv  = new double[numKnots];
//
//	/* read int knots one by one */
//	for(i = 0; i < numKnots; i++)
//	{
//		file.getline(temp, 100);
//		dbl=strtod(temp, &pEnd);
//		dKv[i] = dbl;
//	}//end for i
//	
//	/* finally close stream */
//	file.close();
//}

void readNURBSDataFromFile(const char *strFileName, vector<CNurbs*> &curves)
{
	
	int i;

	/* open stream for reading */
	ifstream file(strFileName);

	if(!(file.is_open()))
	{
		cerr<<"Could not open file "<<strFileName<<endl;
		return;
	}

	char temp[100];
	char *pEnd;
	double dbl;

	VECTOR2 *vCp;
	Real *dKv;
	Real *weights;

	/* read number of objects */
	file.getline(temp, 100);
	int numObjects = atoi(temp);

		
	for(int j = 0; j < numObjects; j++)
	{
		/* read number of points */
		file.getline(temp, 100);
		int numPoints = atoi(temp);


		/* allocate memory for array of control points */
		vCp = new VECTOR2[numPoints]; 
		//cout<<"reading points"<<endl;
		/* read in control points one by one */
		for(i = 0; i < numPoints; i++)
		{
			file.getline(temp, 100);
			dbl=strtod(temp, &pEnd);
			vCp[i].x = dbl;
			file.getline(temp, 100);
			dbl=strtod(temp, &pEnd);
			vCp[i].y = dbl;
			//cout<<vCp[i];
			
		}//end for i

		/* read number of weights */
		file.getline(temp, 100);
		numPoints = atoi(temp);

		/* allocate memory for the weights */
		weights = new Real[numPoints];

		/* read int weights one by one */
		for(i = 0; i < numPoints; i++)
		{
			file.getline(temp, 100);
			dbl=strtod(temp, &pEnd);
			weights[i] = dbl;
		}//end for i


		/* read number of knots */
		file.getline(temp, 100);
		int numKnots = atoi(temp);

		/* allocate memory for knot vector */
		dKv  = new Real[numKnots];

		/* read int knots one by one */
		for(i = 0; i < numKnots; i++)
		{
			file.getline(temp, 100);
			dbl=strtod(temp, &pEnd);
			dKv[i] = dbl;
		}//end for i

		curves.push_back(new CNurbs(3, numPoints, dKv, weights, vCp, 3));

		delete[] dKv;
		delete[] weights;
		delete[] vCp;
	}//end for j
	/* finally close stream */
	file.close();

}

//void readPointsFromFile(const char *strFileName, vector<VECTOR2> &testPoints)
//{
//
//	/* open stream for reading */
//	ifstream file(strFileName);
//
//	if(!(file.is_open()))
//	{
//		cerr<<"Could not open file "<<strFileName<<endl;
//		return;
//	}
//
//	char temp[100];
//	char *pEnd;
//	double dbl;
//	cout <<"Reading points..."<<endl;	
//	/* read in control points one by one */
//	while(!file.eof())
//	{
//		VECTOR2 tmp;
//		file.getline(temp, 100);
//		dbl=strtod(temp, &pEnd);
//		tmp.x = dbl;
//		file.getline(temp, 100);
//		dbl=strtod(temp, &pEnd);
//		tmp.y = dbl;
//		testPoints.push_back(tmp);
//		
//	}//end while
//	cout <<"read "<<testPoints.size()<<" points..."<<endl;
//
//}

void readObjectVertices(const char *sFileName, vector<VECTOR2> &vPoints)
{

	ifstream fin(sFileName);

	int parameters[4];

	int i;

	for(i = 0; i < 4; i++)
	{
		fin >> parameters[i];
	}//end for

	for(i = 0; i < parameters[0]; i++)
	{ 
		int num;
		double dX,dY;
		fin >> num;
		fin >> dX;
		fin >> dY;

		VECTOR2 vec(dX,dY);
		vPoints.push_back(vec);		
		
	}//for
	
	fin.close();

	cout << "Number of object vertices: "<<parameters[0]<<endl;
	
}//end readObjectVertices

void readVoronoiVertices(const char *sFileName, vector<VECTOR2> &vAllVertices)
{

	ifstream fin(sFileName);

	int parameters[4];

	int i;
	
	
	for(i = 0; i < 4; i++)
	{
		fin >> parameters[i];
	}//end for

	for(i = 0; i < parameters[0]; i++)
	{ 
		int num;
		double dX,dY;
		fin >> num;
		fin >> dX;
		fin >> dY;

		
        VECTOR2 vec(dX,dY);
		
		vAllVertices.push_back(vec);		
		
	}//end for
	
	fin.close();

}//end readVoronoiVertices

void readVoronoiEdges(const char* sFileName)
{

	ifstream fin(sFileName);

	if(!fin.is_open())
	{
		cout<<"error opening file "<<sFileName<<endl;
		exit(0);
	}

	int parameters[2];

	int i;
		
	for(i = 0; i < 2; i++)
	{
		fin >> parameters[i];
	}//end for

	for(i = 0; i < parameters[0]; i++)
	{ 
		int num;
		int iStart;
		int iEnd;
		double dX,dY;

		fin >> num;
		fin >> iStart;
		fin >> iEnd;

		if(iEnd == -1)
		{
			fin >> dX;
			fin >> dY;
			continue;
		}

		iStart -=1;
		iEnd   -=1;
		
	}//for
	
	fin.close();

}//end readVoronoiEdges