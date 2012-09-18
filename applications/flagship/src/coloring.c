#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*******************************************************************************
 * The macro FNAME converts a C function name to the form that is
 * required when it is called by a FORTRAN program. On many UNIX
 * systems, subprograms called from FORTRAN have an implied underscore
 * character at the ends of their names. This macro takes care of this
 * operating system quirk.
 *******************************************************************************/
#ifdef VMS
#define FNAME(name)	name
#else
#ifdef __APPLE__
#define FNAME(name)	name
#else
#ifdef __STDC__
#define FNAME(name)	name##_
#else
#define FNAME(name)	name/**/_
#endif
#endif
#endif

#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )

#ifndef STOPWATCH_GUARD
#define STOPWATCH_GUARD

void startTimer (void);		// start the stopwatch

double stopTimer (void);	// stop time measurement and return elapsed time in seconds

#endif  // define STOPWATCH_GUARD

#include <sys/time.h>
#include <stdlib.h>

static struct timeval timer_start; // global variable to store the starting time of the stopwatch


void startTimer (void)		// start the stopwatch
{
   gettimeofday(&timer_start, NULL);
}


double stopTimer (void)		// stop time measurement and return duration in seconds
{
   struct timeval end;
   double elapsed;
   gettimeofday(&end, NULL);
   elapsed = (end.tv_sec-timer_start.tv_sec) + (end.tv_usec -
timer_start.tv_usec)/1000000.0;   
return elapsed;
}




typedef struct {				
  int color;
  int node[2];
  int data[4];
} Edge;

//Transformierte struct
typedef struct {					
  int degree;
	int index;
	int neuergrad;
	int *nachbar;
	int zwischenspeicher;
	int pfad;
} Vertice;

//Adjazenzstruct
typedef struct {	
	struct node *kanteknoten;				
  //int adjazenzanz;
	//int *adjazenz;
	int miscol;
	//int mark;
} Knoten;


typedef struct {	
	struct node *liste;				
} Adjazenz;

typedef struct {	
	struct node *pointer;				
} Pointer;

//struct zum ueberpruefen
typedef struct {				
  int pruefdegree;
	int *pruefnachbarn;
} Pruefer;



struct node {
	struct node *left;
	struct node *right;
	int Kante;
	int Knoten;
};


struct node* new_list(Knoten *knotenliste, int dataKnoten1,int dataKnoten2, int dataKante)	
{
	struct node *new = (struct node*) malloc(sizeof(struct node));
	
	new->Knoten	= dataKnoten2;
	new->Kante	= dataKante;
	new->right = NULL;
	new->left  = NULL;
	return new;
}

struct node* insert(Knoten *knotenliste, int dataKnoten1,int dataKnoten2, int dataKante)								
{										
	struct node *new = (struct node *) malloc(sizeof(struct node));
	
	new->Knoten	= dataKnoten2;
	new->Kante	= dataKante;
	new->right  = knotenliste[dataKnoten1].kanteknoten ;
	new->left   = NULL;
	knotenliste[dataKnoten1].kanteknoten->left  = new;
	return new;
}


void delete(Knoten *gradliste, int Grad ,struct node* Kante)									
{
	//mind. 2 El.  und ganz hinten
	if( Kante->right == NULL  && Kante->left != NULL)		
		{	//printf("Fall1 \n" );
			Kante->left->right=NULL;
			free(Kante);
		}	
	//mind. 2 El und ganz vorne
	else if(Kante->left == NULL  && Kante->right != NULL) 
		{	//printf("Fall2 \n" );
			Kante->right->left=NULL;
			gradliste[Grad].kanteknoten=Kante->right;
			free(Kante);
		}
	//nur ein Element
	else if(Kante->right == NULL && Kante->left == NULL) 
		{	//printf("Fall3 \n" );
			gradliste[Grad].kanteknoten=NULL;
			free(Kante);
		}
	// mind. 2 El. und in der Mitte
	else if(Kante->right != NULL && Kante->left != NULL) 
		{	//printf("Fall4 \n" );
			Kante->left->right=Kante->right;
			Kante->right->left=Kante->left;	
			free(Kante);
		}  
}

void umsetzen(Knoten *gradliste, Pointer *pointerliste, int Grad, int Kante)									
{
		//mind. 2 El.  und ganz hinten
	if( pointerliste[Kante].pointer->right == NULL  && pointerliste[Kante].pointer->left != NULL)		
		{	
			if( gradliste[Grad-1].kanteknoten != NULL)		
				{
					//printf("Fall1.1 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->left->right=NULL;
					//in anderer liste
					pointerliste[Kante].pointer->right=gradliste[Grad-1].kanteknoten;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
					pointerliste[Kante].pointer->left=NULL;
				}
			else		
				{
					//printf("Fall1.2 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->left->right=NULL;
					//in anderer liste
					pointerliste[Kante].pointer->right=NULL;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
					pointerliste[Kante].pointer->left=NULL;
				}
			
		}	
	//mind. 2 El und ganz vorne
	else if(pointerliste[Kante].pointer->left == NULL  && pointerliste[Kante].pointer->right != NULL) 
		{	
			if( gradliste[Grad-1].kanteknoten != NULL)		
				{
					//printf("Fall2.1 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->right->left=NULL;
					gradliste[Grad].kanteknoten = pointerliste[Kante].pointer->right;
					//in anderer liste
					pointerliste[Kante].pointer->right=gradliste[Grad-1].kanteknoten;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
				}
			else		
				{
					//printf("Fall2.2 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->right->left=NULL;
					gradliste[Grad].kanteknoten = pointerliste[Kante].pointer->right;
					//in anderer liste
					pointerliste[Kante].pointer->right=NULL;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
				}
		}
	//nur ein Element
	else if(pointerliste[Kante].pointer->right == NULL && pointerliste[Kante].pointer->left == NULL) 
		{	
			if( gradliste[Grad-1].kanteknoten != NULL)		
				{
					//printf("Fall3.1 \n" );
					//in eigener liste
					gradliste[Grad].kanteknoten=NULL;
					//in anderer liste
					pointerliste[Kante].pointer->right=gradliste[Grad-1].kanteknoten;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;

				}
			else		
				{
					//printf("Fall3.2 \n" );
					//in eigener liste
					gradliste[Grad].kanteknoten=NULL;
					//in anderer liste
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
				}
			
		}
	// mind. 2 El. und in der Mitte
	else if(pointerliste[Kante].pointer->right != NULL && pointerliste[Kante].pointer->left != NULL) 
		{	
			if( gradliste[Grad-1].kanteknoten != NULL)		
				{
					//printf("Fall4.1 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->left->right=pointerliste[Kante].pointer->right;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer->left;
					//in anderer liste
					pointerliste[Kante].pointer->right=gradliste[Grad-1].kanteknoten;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
					pointerliste[Kante].pointer->left=NULL;
				}
			else		
				{
					//printf("Fall4.2 \n" );
					//in eigener liste
					pointerliste[Kante].pointer->left->right=pointerliste[Kante].pointer->right;
					pointerliste[Kante].pointer->right->left=pointerliste[Kante].pointer->left;
					//in anderer liste
					pointerliste[Kante].pointer->right=NULL;
					gradliste[Grad-1].kanteknoten=pointerliste[Kante].pointer;
					pointerliste[Kante].pointer->left=NULL;
				} 	
		}  

}

/*
 * Wrapper
 */

int edgeColoring (int neq, int nedge, Edge *edgelist)
{ 

	srand (time(NULL));
	int j,k,q;
	int Farbanz;

																									
	int i,c,M;
	int w,z,t,r, temp, tester1, nachbarzaehler,Farbe,h, Anzahl, merker,Erhoeher,m,s, index, falsch,w2,Streicherzaehler,outgrad, nvertice,n1,n2,g,marker, color, zaehler,xi,v,ws,gesamt, miscolzaehler;
	
	Vertice *hilfspointer= (Vertice*) malloc(sizeof(Vertice));
	double zeit;


//Waehle Optioenen fuer Programmdurchlauf

// 1 Kanteninfos // 2 Knotennachbarn und Anz. // 3 Adjazenzliste	// 4 nichts //
int preprocessing=4;

// 1 Largest First // 2 Smallest Last // nichts //
int Sortierung=3;	

// 1 Greedy // 2 Least Used // 3 Random // 4 Match // 5 Match Sl // 6 NTL // 7 nichts //
int Algorithmus=5;

//	1 Kanteninfos	// 2 Knoteninfos	// 3 Kantenanz. in Farben // 4 Zulaessigkeit // 5 nichts
int postprocessing=6;




int Farbpool;
int max=20;									//max Anzahl der Farben
int Farbanzahl[max];
int zulaessig[max];
int Random[max];




Farbpool=4;
 

for (j=0; j<max; ++j)					//Arrays mit 0 initialisieren
	{
		Farbanzahl[j]=0;
		zulaessig[j]=0;
		Random[j]=0;
	}

//Pruefstruktur
Pruefer *prueferlist = (Pruefer*) malloc(nedge*sizeof(Pruefer));

//Knotenstruktur1			
Vertice *verticelist = (Vertice*) malloc(nedge*sizeof(Vertice)); 		
	
printf("\nWerte initialisieren..\n" );
	for (j=0; j<nedge; ++j) 																
		{
    verticelist[j].degree = 0;
		verticelist[j].neuergrad = 0;
		verticelist[j].zwischenspeicher = 0;
		verticelist[j].index = j;
		prueferlist[j].pruefdegree = 0;
		}


	

printf("Anzahl der Knoten berechnen..\n" );  
nvertice=-1;
for (j=0; j<nedge; ++j)																			
	{
		n1=edgelist[j].node[0];
		n2=edgelist[j].node[1];
		if ( nvertice<n1 || nvertice<n2 )			
			{									
					nvertice = max(n1,n2);				
			}
	}
nvertice++;	//weil der Knoten 0 auch ex.




//Knotenstruktur2
Knoten *knotenliste=(Knoten*) malloc(nvertice*sizeof(Knoten)); 	
//outgrad=-1; 


///NTL preprocessing

int *Anfangspos = (int*) malloc( nvertice*sizeof(int) );
int *aktpos = (int*) malloc( nvertice*sizeof(int) );
int *hilfsarray = (int*) malloc( nvertice*sizeof(int) );


typedef struct {				
  int Kante;
	int Knoten;
} Tabelle;

Tabelle *Kantenarray=(Tabelle*) malloc((2*nedge)*sizeof(Tabelle));

int *Pfadarray = (int*) malloc( nedge*sizeof(int) );

//miscol, hilfsarray und Pfadarray initialisieren
for (t=0; t < nvertice  ; ++t)		
	{	
		knotenliste[t].miscol=0;
		hilfsarray[t]=0;
		Pfadarray[ t ] = -1;
	}

//Tabelle erstellen
for (j=0; j<nedge; ++j)											
  {
		hilfsarray[ edgelist[j].node[0] ]++;
		hilfsarray[ edgelist[j].node[1] ]++;
	}

Anfangspos[0]=1;
aktpos[0]=1;	
for (t=1; t < nvertice  ; ++t)  
	{	
		Anfangspos[t]=  Anfangspos[t-1]+hilfsarray[t-1]  ;
		aktpos[t]= Anfangspos[t]  ;		
	}


outgrad=0;
for (t=0; t < nvertice  ; ++t)  
	{	
		outgrad = max(outgrad, hilfsarray[t]);
	}

//nur zum testen da
printf(" outgrad ist %d \n\n" , outgrad);

int colors;
colors= outgrad+1;

int W[colors];
int usedmiscol[colors];
int Fan[colors];
int Endknoten[colors];
int miscolEndknoten[colors];

//Pointer der Arrayliste initialisieren
for (j=0; j<nvertice; ++j)																			
	{
		knotenliste[j].kanteknoten=NULL;
	}


int Knoten1,Knoten2;

if( Algorithmus !=6 || postprocessing == 4 || postprocessing == 2 || Sortierung == 2)													
	{

startTimer ();	

//degree aus hilfsarray bestimmen
for (j=0; j<nedge; ++j)																			
	{
		verticelist[j].degree=hilfsarray[ edgelist[j].node[0] ]+hilfsarray[ edgelist[j].node[1] ]-2;
		verticelist[j].neuergrad=verticelist[j].degree;
		prueferlist[j].pruefdegree=verticelist[j].degree;	
	}	

Knoten *adjazenzliste=(Knoten*) malloc(nvertice*sizeof(Knoten)); 

//adjazenzlistnepointer initialisieren
for (j=0; j<nvertice; ++j)																			
	{
		adjazenzliste[j].kanteknoten=NULL;
	}

//adjazenliste bestimmen
for (j=0; j<nedge; ++j)
    {

			Knoten1=edgelist[j].node[0];
			Knoten2=edgelist[j].node[1];
			
			//fuege die kante j in die adjazenzliste ein
			if ( adjazenzliste[Knoten1].kanteknoten == NULL )  
				{
					adjazenzliste[Knoten1].kanteknoten=new_list(adjazenzliste,Knoten1,Knoten2,j);
				}
			else
				{
					adjazenzliste[Knoten1].kanteknoten=insert(adjazenzliste,Knoten1,Knoten2,j);
				}

			if ( adjazenzliste[Knoten2].kanteknoten == NULL )  
				{
					adjazenzliste[Knoten2].kanteknoten=new_list(adjazenzliste,Knoten2,Knoten1,j);
				}
			else
				{
					adjazenzliste[Knoten2].kanteknoten=insert(adjazenzliste,Knoten2,Knoten1,j);
				}
	}


printf("Speicher fuer Nachbarn allokieren..\n" );	
for (j=0; j<nedge; ++j)												
	{
		verticelist[j].nachbar = (int*) malloc(  (verticelist[j].degree) *sizeof(int));
		prueferlist[j].pruefnachbarn = (int*) malloc((prueferlist[j].pruefdegree) *sizeof(int));
	}

//nachbarn aus der adjazenzliste entnehmen
struct node *writer;
for (j=0; j<nedge; ++j)
  {
		g=-1;

		Knoten1=edgelist[j].node[0];
		Knoten2=edgelist[j].node[1];

		writer= adjazenzliste[Knoten1].kanteknoten;
		//printf("j ist %d" , j);
		while ( writer != NULL   )
			{	//printf("Kante %d // " ,writer->Kante);	
				
				if( writer->Kante != j)
					{	
						g++;
						verticelist[j].nachbar[g]= writer->Kante;
						prueferlist[j].pruefnachbarn[g]= writer->Kante;	
								
					}
				writer=writer->right;
				
			}
		
		writer= adjazenzliste[Knoten2].kanteknoten;
		while ( writer!= NULL   )
			{	//printf("Kante %d // " ,writer->Kante);	
				
				if(writer->Kante != j)
					{	
						g++;
						verticelist[j].nachbar[g]= writer->Kante;
						prueferlist[j].pruefnachbarn[g]= writer->Kante;	
								
					}
				writer=writer->right;
			}
		//printf("\n\n\n" );
	}

zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);
}



startTimer ();

//pointer der arrayliste initialisieren
Pointer *pointerliste = (Pointer*) malloc(nedge*sizeof(Pointer)); 
for (j=0; j<nedge; ++j)
  {
		pointerliste[j].pointer=NULL;
	}


//minimalen und maximalen grad bestimmen
int maxgrad;
int mingrad;
mingrad=999;
maxgrad=0;
for (j=0; j<nedge; ++j)
  {
		mingrad=min(verticelist[j].degree, mingrad);
		maxgrad=max(verticelist[j].degree, maxgrad);
	}
printf(" mingrad ist %d\n" , mingrad);
printf(" maxgrad ist %d\n" , maxgrad);


//speicher fuer gradarray allokieren
Knoten *gradliste=(Knoten*) malloc((maxgrad+1)*sizeof(Knoten));


//gradlistenpointer initialisieren
for (j=0; j<(maxgrad+1); ++j)																			
	{
		gradliste[j].kanteknoten=NULL;
	}

int gradspeicher;
//gradliste bestimmen
for (j=0; j<nedge; ++j)
    {
	
			gradspeicher=verticelist[j].degree;

			//fuege die kante j in die gradliste ein
			if ( gradliste[gradspeicher].kanteknoten == NULL )  
				{
					gradliste[gradspeicher].kanteknoten=new_list(gradliste,1,gradspeicher,j);
					pointerliste[j].pointer= gradliste[gradspeicher].kanteknoten;
				}
			else
				{
					gradliste[gradspeicher].kanteknoten=insert(gradliste,gradspeicher,gradspeicher,j);
					pointerliste[j].pointer= gradliste[gradspeicher].kanteknoten;
				}
	}

zeit= stopTimer ();
printf(" Zeit von gradliste ist %f \n\n" , zeit);

/*
//nur zum testen da
for (j=0; j<nedge; ++j)
  {
		printf(" / %d /" , pointerliste[j].pointer->Kante);
	}
*/

//positionsarray fuer sortierungen allokieren und initialisieren
int *positionsarray = (int*) malloc( nedge*sizeof(int) );
for (j=0; j<nedge; ++j)
  {
		positionsarray[j]=j;
	}

//streicherarray fuer smallest last sortierung allokieren und initialisieren
int *streicherarray = (int*) malloc( nedge*sizeof(int) );
for (j=0; j<nedge; ++j)
  {
		streicherarray[j]=0;
	}



/*
if( (Algorithmus !=4 && Algorithmus !=6) ||postprocessing == 4 || postprocessing == 2)													
	{
startTimer ();

printf("Anzahl der Nachbarn bestimmen..\n" );	
	for (j=0; j<nedge; ++j)									
	{
		verticelist[j].degree=0;		
		for (k=0; k<nedge; ++k)
 		{
			if (  j!=k && ( edgelist[j].node[0]==edgelist[k].node[0] || edgelist[j].node[0]==edgelist[k].node[1] || edgelist[j].node[1]==edgelist[k].node[0] || edgelist[j].node[1]==edgelist[k].node[1])   )
			{
				verticelist[j].degree++;
				verticelist[j].neuergrad++; 
				prueferlist[j].pruefdegree++;
			}
		}
	}


printf("Speicher fuer Nachbarn allokieren..\n" );	
for (j=0; j<nedge; ++j)												
	{
		verticelist[j].nachbar = (int*) malloc(  (verticelist[j].degree) *sizeof(int));
		prueferlist[j].pruefnachbarn = (int*) malloc((prueferlist[j].pruefdegree) *sizeof(int));
	}

printf("Nachbarn bestimmen..\n\n" );
g=-1;
	for (j=0; j<nedge; ++j)															
	{
		for (k=0; k<nedge; ++k)
 		{
			if (j!=k && ( edgelist[j].node[0]==edgelist[k].node[0] || edgelist[j].node[0]==edgelist[k].node[1] || edgelist[j].node[1]==edgelist[k].node[0] ||  edgelist[j].node[1]==edgelist[k].node[1]))
			{
			g++;
			verticelist[j].nachbar[g]=k;
			prueferlist[j].pruefnachbarn[g]=k;
			}
		}
		g=-1;
	}

zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

}
*/

switch(preprocessing)
{
case 1:  //Kanteninfos ausgeben//

	for (j=0; j<nedge; ++j)
    {
		printf("Kante %d hat Knoten %d %d und Farbe %d\n", j, edgelist[j].node[0], edgelist[j].node[1], edgelist[j].color);
 		}

break;

case 2: //Knotennachbarn ausgeben//

for (j=0; j<nedge; ++j)
    {
			printf("Knoten %d hat die Nachbarn:", j);
			for (nachbarzaehler=0; nachbarzaehler<verticelist[j].degree; ++nachbarzaehler)
			{
				printf(" %d ", verticelist[j].nachbar[nachbarzaehler]);
			}
		printf("\n");
		}	printf("\n");

				//Knotengrad ausgeben//

	for (j=0; j<nedge; ++j)
    {
		printf("Knoten %d hat Grad %d\n", j, verticelist[j].degree);
		}	printf("\n");

break;

case 3: //Adjazenzliste ausgeben//

printf("\n");
struct node *look;
for (j=0; j<(maxgrad+1); ++j)
	{	
		look= gradliste[j].kanteknoten;
		printf("Knoten %d: ",j );
		while ( look!= NULL   )
		{	printf("Kante %d // " ,look->Kante);	
			look=look->right;
		}printf("\n");
	}


break; 
case 4: //nichts//		
break;
}

 


switch(Sortierung)
{
case 1:				/////Largest First////
printf("/////////////LARGEST FIRST//////////////\n\n");
printf("Anzahl der Kanten: %d\n\n", nedge);


startTimer ();

struct node *look2;
struct node *tempKante;
int aktKante;
//printf(" maxgrad ist %d \n" , maxgrad);

j=0;



	while (j<nedge)
		{
			look2=gradliste[maxgrad].kanteknoten;
			if( look2 != NULL  )
				{
					//startTimer ();
					aktKante=look2->Kante;
					tempKante=look2;		
					//while ( look2 != NULL   )
						//{	//printf("%d // " ,look2->Kante);	
							
							//if( aktKante > look2->Kante )
							//	{
								//	aktKante=look2->Kante;
								//	tempKante=look2;
								//}
							//look2=look2->right;
						//}
					
				//printf("\n\n Loesche %d \n" ,aktKante);
				//printf("zeige auf %d \n" ,tempKante->Kante);

					
					delete(gradliste, maxgrad ,tempKante);
					
					if( j!= positionsarray[aktKante])
						{		
			//printf("\n Tausche! j ist %d \n\n", j );
	//printf("Sind bei Position %d, kante %d muss dahin\n", j, aktKante );
	//printf("Position von Kante %d ist %d\n\n ", aktKante,positionsarray[aktKante]);
	//printf("Auf der Position %d befindet sich die Kante %d\n", j, verticelist[j].index );



		//printf("\n positionsarray vor dem tausch: ");
		//for (h=0; h<nedge; ++h)
    //{
		//	printf(" %d ", positionsarray[h]);
		//}printf("\n\n" );
		

							*hilfspointer =  verticelist[j];
							verticelist[j]	=	verticelist[positionsarray[aktKante]];
							verticelist[positionsarray[aktKante]] =	*hilfspointer	 ;
							
positionsarray[verticelist[positionsarray[aktKante]].index]=	positionsarray[aktKante];
							positionsarray[aktKante]=	j;

				
			
		//for (h=0; h<nedge; ++h)
    //{
		//printf("Knoten %d hat Grad %d und alten Index %d\n", h, verticelist[h].degree, verticelist[h].index);
		//}

		//printf(" positionsarray nach dem tausch: ");
		//for (h=0; h<nedge; ++h)
    //{
		//	printf(" %d ", positionsarray[h]);
		//}
		
						}
					//else
						//{	//printf("\n kein Tausch noetig! j ist %d \n" , j);
						//}
					j++;
				}
			else
				{	//printf("\n\n" );
					maxgrad--;
					//printf("Liste ist leer!!\n" );
				}
			
		}



/*
		for (j=0; j<nedge; ++j)
			{
				z= verticelist[j].degree;
				index=verticelist[j].index;
				m=0;
				for (k=j+1; k<nedge; ++k)
					{
					if ( verticelist[k].degree > z || ( verticelist[k].degree == z && 
							index > verticelist[k].index))
						{
						z= verticelist[k].degree;
						index=verticelist[k].index;
						m++;
						w=k;
						}
					}
				if(m>0)
					{		
						*hilfspointer =  verticelist[j];
						verticelist[j]	=	verticelist[w];
						verticelist[w] =	*hilfspointer	 ;		
					}
		}


*/

zeit= stopTimer ();
printf(" \n\nZeit ist %f \n\n" , zeit);

break;

case 2: //////Smallest Last/////
printf("/////////////SMALLEST LAST//////////////\n\n");
printf("Anzahl der Kanten: %d\n\n", nedge);


startTimer ();



///*
struct node *look3;
struct node *tempKante2;
int aktKante2;
//printf(" maxgrad ist %d \n" , maxgrad);
int nachbar;
j=nedge-1;

	while (j>0)
		{
			look3=gradliste[mingrad].kanteknoten;
			if( look3 != NULL  )
				{					
					aktKante2=look3->Kante;
					tempKante2=look3;		
						
					printf("Delete %d \n" ,aktKante2);			
					delete(gradliste, mingrad ,tempKante2);
					streicherarray[aktKante2]=1;

					
					////nur zum testen da
					//printf("streicherarray: " );
					//for (k=0; k<nedge; ++k)
 					// {
					//	printf(" %d " ,streicherarray[k]);
					//}printf("\n\n" );
					

					//bis jetzt nur zum testen da
					//printf("Die Nachbarn von %d sind: " ,aktKante2);
					for (t=0; t<verticelist[positionsarray[aktKante2]].degree; ++t)
						{			//printf(" %d " , verticelist[positionsarray[aktKante2]].nachbar[t]);
									nachbar= verticelist[positionsarray[aktKante2]].nachbar[t];
							if( streicherarray[ nachbar ] == 0  )
								{	//printf("Der Grad von der Kante ist %d: " ,pointerliste[nachbar].pointer->Knoten);
									umsetzen(gradliste, pointerliste,pointerliste[nachbar].pointer->Knoten , nachbar);	
									pointerliste[nachbar].pointer->Knoten--;
									mingrad= min(mingrad, (pointerliste[nachbar].pointer->Knoten));
							//printf("\n////mingrad ist %d://// \n" ,mingrad);

										
									//struct node *look10;
									//for (h=0; h<(maxgrad+1); ++h)
									//	{	
									//		look10= gradliste[h].kanteknoten;
									//		printf("Knoten %d: ",h );
									//		while ( look10 != NULL   )
									//			{	printf("Kante %d // " ,look10->Kante);	
									//				look10 =look10->right;
									//			}printf("\n");
									//	}printf("\n\n\n\n");
									//	
								}
							
								
						}	//printf("\n\n" );



					if( j!= positionsarray[aktKante2])
						{		
							*hilfspointer =  verticelist[j];
							verticelist[j]	=	verticelist[positionsarray[aktKante2]];
							verticelist[positionsarray[aktKante2]] =	*hilfspointer	 ;
							
							positionsarray[verticelist[positionsarray[aktKante2]].index]=	positionsarray[aktKante2];
							positionsarray[aktKante2]=	j;
						}
					else
						{	//printf("\n kein Tausch noetig! j ist %d \n" , j);
							
						}
					j--;
				}
			else
				{	//printf("\n\n" );
					mingrad++;
					//printf("Liste ist leer!!\n" );
				}
			
		}

//*/

/*
w=0;
index=0;
	for (j=nedge-1; j>0; j--)
		{
			z= verticelist[j].neuergrad; 
			index=verticelist[j].index;  		
			m=0;
			for (k=j-1; k>-1; k--)
				{
					if (verticelist[k].neuergrad < z  || ( verticelist[k].neuergrad == z && index > 							verticelist[k].index))		
						{//printf("Wert1: %d Wert2: %d Wert3: %d Wert4: %d\n", verticelist[k].neuergrad ,z , index, verticelist[k].index);
						z= verticelist[k].neuergrad;			
						index=verticelist[k].index;
						w=k;
						m++;
						}
				}
			if(m>0)		//geandert
				{	//printf("Tausch: %d mit %d \n", j,w);
					*hilfspointer =  verticelist[j];
					verticelist[j]	=	verticelist[w];
					verticelist[w] =	*hilfspointer	 ;
				}
				else ; //printf("%d kein Tausch!\n", j);
					for (t=j-1; t>-1; t--)
					{
						for (r=0; r<verticelist[t].neuergrad; r++)		
						{
							if (verticelist[t].nachbar[r]==verticelist[j].index)
							{//printf("Streiche in %d die %d !!!\n", t ,verticelist[t].nachbar[r]);
								temp=verticelist[t].nachbar[r];
								verticelist[t].nachbar[r]=verticelist[t].nachbar[verticelist[t].neuergrad-1]; 
								verticelist[t].nachbar[verticelist[t].neuergrad-1]= temp;	 										
								verticelist[t].neuergrad--;																										
							}
						}
					}
						//printf("\n");
						//printf("///////Nach dem %d/ten Tausch/////////\n\n", nedge-j);
						//for (tester1=0; tester1<nedge; ++tester1)
    				//	{
			//printf("Knoten %d mit altem index: %d und Grad: %d  \n", tester1,verticelist[tester1].index, verticelist[tester1].neuergrad);
						//for (nachbarzaehler=0; nachbarzaehler<11; ++nachbarzaehler)
							//{
							//	printf(" %d ", verticelist[tester1].nachbar[nachbarzaehler]);
							//}
							//printf("\n");
							//}
						//printf("\n\n");	 	
					}

*/

zeit= stopTimer ();
printf(" \n\nZeit ist %f \n\n" , zeit);

break;

case 3:	// nichts //
break;
}



switch(Algorithmus)
{ 

case 1: 	//Greedy//
printf("/////////////GREEDY//////////////\n\n");
printf("Anzahl der Kanten: %d\n\n", nedge);

startTimer ();

Farbanz=0;
Farbe=0;
	
for (j=0; j<nedge; ++j)								
	{	

//if(j==5000 || j==10000 || j==15000 || j==20000 || j==25000 || j==30000 || j==35000|| j==40000 || j==45000 || j==50000 || j==55000 || j==60000 || j==65000)
	//{	printf(" j ist %d\n",j);
	//}


		//solange ungefaerbt suche Farbe 
		while (edgelist[	verticelist[j].index	].color==0)  
			{
				Farbe++;
				g=0;
				//Nachbarn durchgehen
				for (k=0; k<verticelist[j].degree; ++k) 					
					{
						//haben Nachbarn aktuelle Farbe? nein, => g++ 
						if (Farbe != edgelist[ verticelist[j].nachbar[k]  ].color)  
							{
								g++; 															
							}
					}
				// =g wenn kein Nachbar akt. Farbe besitzt
				if (g==verticelist[j].degree)							
					{
							edgelist[ verticelist[j].index ].color=Farbe;	
							Farbanzahl[Farbe-1]++;			
							//printf("Farbe von Kante %d: %d\n", verticelist[j].index, Farbe);
					}
			}
		if (Farbe>Farbanz)				
			{
				Farbanz=Farbe;
			}
	Farbe=0;
	}	//printf("Farbanzahl ist %d\n", Farbanz);


zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

break;

case 2: ///Greedy Least Used////

printf("/////////GREEDY LEAST USED//////////////\n\n");
printf("Anzahl der Knoten: %d\n\n", nedge);

startTimer ();

g=0;
merker=0;
	
	for (j=0; j<nedge; ++j)
	{
		//Schleife ueber die versch. Farben
		for (k=0; k<Farbpool; ++k)
			{ g=0;
				//schleife ueber die nachbarn
				for (t=0; t<verticelist[j].degree; ++t)
					{
						//haben Nachbarn aktuelle Farbe? nein, => g++
						if ( k+1 != edgelist[ verticelist[j].nachbar[t]   ].color)
							{
								g++; 							
							}
					}		
				//falls Kante gefaerbt und g=Knotengrand und akt. Anzahl kleiner als gemerkte faerbe um und merke dir neue Position
				if ( edgelist[  verticelist[j].index    ].color!=0 && g==verticelist[j].degree && Farbanzahl[k]<Farbanzahl[merker])   
					{	
						merker=k;
						edgelist[ verticelist[j].index ].color=k+1;			
					}
				//falls Kante ungefaerbt und g=Knotengrad faerbe Kante und merke dir Position
				else if ( edgelist[   verticelist[j].index  ].color==0 && g==verticelist[j].degree)
					{			
						merker=k; 		
						edgelist[  verticelist[j].index   ].color=k+1;												
					}
			}
		//wenn Kante gefaerbt wurde erhoehe Farbanzahl an passender Stelle
		if( edgelist[   verticelist[j].index   ].color != 0)				
			{
				Erhoeher=edgelist[  verticelist[j].index    ].color-1;				
				Farbanzahl[Erhoeher]++;
			}
		//wenn Kante ungefaerbt ist erhoehe den Farbpool, waehle fuer Kante die neue Farbe
		else
			{
	//printf("Zu wenig farben! Ich erhoehe den Farbpool auf %d in Schritt %d!\n\n",Farbpool+1,j);	
				Farbpool++;
				edgelist[   verticelist[j].index   ].color=Farbpool;			
				Farbanzahl[Farbpool-1]++;
				//printf("Ich faerbe mit %d !!!!!!!!!!!!\n\n" ,Farbpool);
			}
		//falls Farbpool max ueberschreitet setze manuell max um und starte das Programm neu!
		if( Farbpool > max)
			{
			printf("Erhoehe max!\n"); 
			break;
			}
}

Farbanz=Farbpool;		//printf("Farbanzahl ist %d\n", Farbanz);


zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

break;

case 3: ///Random2///
printf("/////////GREEDY RANDOM//////////////\n\n");
printf("Anzahl der Knoten: %d\n\n", nedge);

startTimer ();

tester1=0;
falsch=0;
h=0;
Farbe=0;
	
for (j=0; j<nedge; ++j)
	{
		h=0;
		//Schleife ueber die Farben
		for (k=0; k<Farbpool; ++k)
			{
				g=0;
				//schleife ueber die Nachbarn
				for (t=0; t<verticelist[j].degree; ++t)
					{
						//haben Nachbarn aktuelle Farbe? nein, => g++
						if ( k+1 != edgelist[verticelist[j].nachbar[t]].color)
						{
							g++; 
						}
					}
				//falls akt. Farbe zulaessig, schreibe sie ins Randomarray
				if (g==verticelist[j].degree && zulaessig[k]==0)
					{
						Random[h]=k+1;
						h++;						
					}
			}

		//falls es eine zulaessige Farbe gibt, waehle zufaellig eine aus
		if (h!= 0)
			{
				tester1=(rand()%h) ;	// kann man auch weglassen und unten einbauen
				Farbe	= Random[ tester1 ];
				//printf("Farbe ist der %d stelle ist %d \n", tester1 ,Farbe);
				edgelist[ verticelist[j].index   ].color= Farbe;
				//printf("Farbe von Knoten %d: %d \n",j, edgelist[ verticelist[j].index].color);
				Farbanzahl[Farbe-1]++; 
			}
		//falls es keine zulaessige gibt, erhoehe Farbpool und faerbe mit neuer Farbe
		else
			{
				Farbpool++;
				edgelist[  verticelist[j].index  ].color=Farbpool; // vorher stand j in der klammer
				Farbanzahl[Farbpool-1]++;
				//printf("Farbe von Knoten %d: %d \n",j, edgelist[ verticelist[j].index].color);	
			}
		//ueberschreibe die Randomarrayeintraege mit Nullen fuer den naechsten Durchlauf!
		for (s=0; s<max; ++s)
			{
				Random[s]=0;
			}
		//printf("\n\n");
	}

zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

Farbanz=Farbpool;		//printf("Farbanzahl ist %d\n", Farbanz);
break;

case 4 :		////Greedy-Match Edgecolor///
printf("////////MATCH-EDGECOLOR///////////\n\n");
printf("Anzahl der Knoten: %d\n\n", nedge);

startTimer ();


i=1;
c=0;
M = nedge;			//Anzahl der Kanten im Restgraphen

//falls Graph noch nicht ganz gefaerbt
while (c<nedge)
{
	for (j=0; j<nedge; ++j)
		{
			//falls Kante nicht gefaerbt, faerbe sie mit aktueller Farbe
			if (edgelist[j].color == -i+1)
				{
					//printf(" Faerbe %d \n" , j);
					edgelist[j].color = i;
					c++;
					Farbanzahl[i-1]++;
					//Streiche die Nachbarn
					for (k=0; k< verticelist[j].degree; ++k)
						{
							//printf(" Nachbar %d \n" , verticelist[j].nachbar[k]);
							if ( edgelist[ verticelist[j].nachbar[k] ].color <1)
							{
								edgelist[ verticelist[j].nachbar[k] ].color = -i ;
							}
						}
				
				}
		}
	//printf(" \n\n");
	//erhoehe Farbe
	i++;
}

/*
i=1;
c=0;
M = nedge;			//Anzahl der Kanten im Restgraphen

//falls Graph noch nicht ganz gefaerbt
while (c<nedge)
{
	//falls Teilgraph noch nicht ganz gefaerbt
	while (M>0)
	{
		//schleife uber die Kanten
		for (j=0; j<nedge; ++j)
		{
			//falls Kante nicht gefaerbt, faerbe sie mit aktueller Farbe
			if (edgelist[j].color==-i+1)
			{
				edgelist[j].color= i;
				Farbanzahl[i-1]++;
				c++;
				M--;
				//schleife ueber die Kanten
				for (k=0; k<nedge; ++k)
				{
					//streiche Nachbarn
					if (  k!=j && (edgelist[k].color<1) && ( edgelist[k].node[0]==edgelist[j].node[0] || edgelist[k].node[0]==edgelist[j].node[1] || edgelist[k].node[1]==edgelist[j].node[0] || edgelist[k].node[1]==edgelist[j].node[1])   )
					{
						edgelist[k].color= -i;
						M--;
					}
				}
			}
		}
		//erhoehe Farbe, und berechne neue Anzahl im teilgraphen
		i++;
		M=nedge-c;
	}
}
*/

zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

 Farbanz=i-1;  //printf("Farbanzahl ist %d\n", Farbanz);
break;



case 5: ///MATCH-EDGECOLOR2///

printf("////////MATCH-EDGECOLOR2!!!!///////////\n\n");


			/*
			printf("Streicherarray:");
			for (r=0; r<nedge; ++r) 
				{
					printf(" %d ", streicherarray[ r ]);
				}printf("\n\n");	
			*/

/*
int nachbar;
j=nedge-1;

	while (j>0)
		{
			look3=gradliste[mingrad].kanteknoten;
			if( look3 != NULL  )
				{					
					aktKante2=look3->Kante;
					tempKante2=look3;		
						
					printf("Delete %d \n" ,aktKante2);			
					delete(gradliste, mingrad ,tempKante2);
					streicherarray[aktKante2]=1;
					//printf("Die Nachbarn von %d sind: " ,aktKante2);
					for (t=0; t<verticelist[positionsarray[aktKante2]].degree; ++t)
						{			//printf(" %d " , verticelist[positionsarray[aktKante2]].nachbar[t]);
									nachbar= verticelist[positionsarray[aktKante2]].nachbar[t];
							if( streicherarray[ nachbar ] == 0  )
								{	//printf("Der Grad von der Kante ist %d: " ,pointerliste[nachbar].pointer->Knoten);
									umsetzen(gradliste, pointerliste,pointerliste[nachbar].pointer->Knoten , nachbar);	
									pointerliste[nachbar].pointer->Knoten--;
									mingrad= min(mingrad, (pointerliste[nachbar].pointer->Knoten));					
								}					
						}	
					if( j!= positionsarray[aktKante2])
						{		
							*hilfspointer =  verticelist[j];
							verticelist[j]	=	verticelist[positionsarray[aktKante2]];
							verticelist[positionsarray[aktKante2]] =	*hilfspointer	 ;
							
							positionsarray[verticelist[positionsarray[aktKante2]].index]=	positionsarray[aktKante2];
							positionsarray[aktKante2]=	j;
						}
					else
						{	//printf("\n kein Tausch noetig! j ist %d \n" , j);						
						}
					j--;
				}
			else
				{	//printf("Liste ist leer!!\n" );
					mingrad++;
				}
			
		}

//*/
startTimer ();

j=0;
i=1;
c=0;
struct node *look4; struct node *tempKante3; int aktKante3; int nachbar2; int nachbar3;

//j ist fuer die gefaerbten Kanten
while (j<nedge)
	{
		// c ist fuer die gefaerbten und gestrichenen Kanten
		
		while (c<nedge)
			{
				//pointer wird auf Kante mit mingrad gesetzt
				look4=gradliste[mingrad].kanteknoten;
				if( look4 != NULL  )
					{					
						aktKante3=look4->Kante;
						tempKante3=look4;	

						//Faerbe Kante u. erhoehe Farbanzahl und aktuell gefae./gestr. Kanten 
						edgelist[ verticelist[aktKante3].index   ].color= i;
						c++;
						j++;
						Farbanzahl[i-1]++;

						//loesche die Kante aus der Gradliste denn sie ist schon gefaerbt
						delete(gradliste, mingrad ,tempKante3);
						streicherarray[aktKante3]=1;							//vllt mit Farbe ersetzen

							//printf("1.Delete %d \n\n" , aktKante3);

									/*
									printf("Druck 1 \n" );
									struct node *look10;
									for (h=0; h<(maxgrad+1); ++h)
										{	
											look10= gradliste[h].kanteknoten;
											printf("Knoten %d: ",h );
											while ( look10 != NULL   )
												{	printf("Kante %d // " ,look10->Kante);	
													look10 =look10->right;
												}printf("\n");
										}	printf("\n\n\n\n");
									*/



						
						for (t=0; t<verticelist[positionsarray[aktKante3]].neuergrad; ++t)
							{	//printf(" %d " , verticelist[positionsarray[aktKante2]].nachbar[t]);
								nachbar2= verticelist[positionsarray[aktKante3]].nachbar[t];
								
								//Reduziere den Grad der Nachbarn fuer den naechsten Durchlauf
								verticelist[nachbar2].degree--;
								if( streicherarray[ nachbar2 ] == 0  )
									{	//printf("Der Grad von der Kante ist %d: " ,pointerliste[nachbar].pointer->Knoten);										
										//loesche Nachbarn 
										delete(gradliste,pointerliste[nachbar2].pointer->Knoten  ,pointerliste[nachbar2].pointer);
										streicherarray[nachbar2]=2;	


											//printf("2.Delete %d \n\n" , nachbar2);
										
										/*
										printf("Druck 2 \n" );
										struct node *look10;
										for (h=0; h<(maxgrad+1); ++h)
											{	
												look10= gradliste[h].kanteknoten;
												printf("Knoten %d: ",h );
												while ( look10 != NULL   )
													{	printf("Kante %d // " ,look10->Kante);	
														look10 =look10->right;
													}printf("\n");
											}	printf("\n\n\n\n");
										*/

										//erhoehe akt. gefae./gestr. Kanten im Graph
										c++;
						
										
										//nachbarn der nachbarn werden gestrichen			
										for (r=0; r<verticelist[positionsarray[nachbar2]].neuergrad; ++r)
											{
												nachbar3= verticelist[positionsarray[nachbar2]].nachbar[r];
												//printf("Nachbar3 ist %d \n\n" , nachbar3);
												if( streicherarray[ nachbar3 ] == 0  )
													{
															//printf("Reduziere %d um 1 \n\n" , nachbar3);
														//nachbarn werden in der Gradliste eins hoeher angesetzt
														umsetzen(gradliste, pointerliste,pointerliste[nachbar3].pointer->Knoten , nachbar3);	
												
														//Grad wird um eins reduziert
														pointerliste[nachbar3].pointer->Knoten--;
										
														//neuer mingrad wird bestimmt
														mingrad= min(mingrad, (pointerliste[nachbar3].pointer->Knoten));
													}
											}			
												
									}					
							}	
					}	
				else
					{	//printf("Liste ist leer!!\n" );
						mingrad++;
					}	
			}

		//printf("\n\n//////////Fertig mit Durchgang//////////////\n\n\n\n" );

		//erhoehe Farbe um 1
		i++;
		c=j;
		/*
		printf("Gradzahlen:\n");
			for (r=0; r<nedge; ++r) 
				{
					printf("Grad von Kante %d ist %d \n",r, verticelist[r].degree );
				}printf("\n\n");
		*/

		/*
		printf("Streicherarray:");
		for (r=0; r<nedge; ++r) 
			{
				printf(" %d ", streicherarray[ r ]);
			}printf("\n\n");
		*/
		
		for (r=0; r<nedge; ++r) 
			{
				if( streicherarray[ r ] == 2  )
					{

						if ( gradliste[ verticelist[r].degree].kanteknoten == NULL )  
							{
								gradliste[ verticelist[r].degree].kanteknoten=new_list(gradliste,1,verticelist[r].degree,r);
								pointerliste[r].pointer= gradliste[ verticelist[r].degree].kanteknoten;
							}
						else
							{
								gradliste[ verticelist[r].degree].kanteknoten=insert(gradliste,verticelist[r].degree,verticelist[r].degree,r);
								pointerliste[r].pointer= gradliste[ verticelist[r].degree].kanteknoten;
							}
						streicherarray[ r ] = 0;

						//printf(" %d",  r  );
					}
			}//printf("\n\n\n\n");

			mingrad= 0;
			/*
			printf("Streicherarray:");
			for (r=0; r<nedge; ++r) 
				{
					printf(" %d ", streicherarray[ r ]);
				}printf("\n\n");	
			*/
		
				/*
				printf("/// Druck 3 ////\n" );
				struct node *look10;
				for (h=0; h<(maxgrad+1); ++h)
					{	
						look10= gradliste[h].kanteknoten;
						printf("Knoten %d: ",h );
						while ( look10 != NULL   )
							{	printf("Kante %d // " ,look10->Kante);	
								look10 =look10->right;
							}printf("\n");
					}	printf("\n\n\n\n");
				*/
	}				

zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);

/*
startTimer ();

i=1;
c=0;
j=0;

//falls Graph noch nicht ganz gefaerbt
while (c<nedge)
{
	//falls Teilgraph noch nicht ganz gefaerbt
	while (j<nedge) 
	{	

		//for (tester1= 0 ; tester1<nedge; ++tester1)    
		//	{
		//		printf("//1111// Kante %d: index %d    Grad %d    Farbe %d   TEMP %d\n", tester1, verticelist[tester1].index, verticelist[tester1].degree, edgelist[ verticelist[tester1].index].color,verticelist[tester1].zwischenspeicher);	
		//	}	printf("\n\n");  


		z= verticelist[j].degree;						
		index=verticelist[j].index;
		w=j;
		//schleife ueber die kanten
		for (k= j ; k<nedge; ++k)     
			{ 
				//falls neue Kante kleineren grad hat tausche(bei gleichheit gucke auf alten index) 
				//merke den Grad und die Position
				if ( ( z > verticelist[k].degree && edgelist[ verticelist[k].index  ].color==-i+1 ) || ( verticelist[k].degree == z  &&  index > verticelist[k].index && edgelist[verticelist[k].index ].color==-i+1)  )						
					{  
						z= verticelist[k].degree;
						index=verticelist[k].index;
						w=k;
					}
			}
		//faerbe gefundene Kante ein
		edgelist[ verticelist[w].index   ].color= i;
		Farbanzahl[i-1]++;
		c++;

if(c==5000 || c==10000 || c==15000 || c==20000 || c==25000 || c==30000 || c==35000|| c==40000 || c==45000 || c==50000 || c==55000 || c==60000 || c==65000  || c==80000 || c==100000 || c==120000 || c==150000 || c==180000 ||  c==200000 ||  c==220000 ||  c==240000 || c==250000 )
	{	printf(" c ist %d\n",c);
	}

		//tausche nun gefaerbte Kante an Stelle die nicht mehr durchgegangen wird
		if(w !=c-1) 
			{ 
				*hilfspointer =  verticelist[c-1];
				verticelist[c-1]	=	verticelist[w];
				verticelist[w] =	*hilfspointer;		
				if( edgelist[ verticelist[w].index].color== -i && w>j)	// kein luecken entstehen!!!
					{		 						
						*hilfspointer =  verticelist[j];
						verticelist[j]	=	verticelist[w];
						verticelist[w] =	*hilfspointer;	
					}
				j++;
			}
		else 
			{
				j++;	
			}

		//Nachbarn werden gestrichen(auch an passende Stelle getauscht)
		for (t=c; t<nedge; ++t) // Streicher
			{ 
				w2=0;
				if (	 (edgelist[verticelist[t].index].color<1) && ( edgelist[ verticelist[t].index].node[0]==edgelist[verticelist[c-1].index].node[0] || edgelist[verticelist[t].index].node[0]==edgelist[verticelist[c-1].index].node[1] || edgelist[verticelist[t].index].node[1]==edgelist[verticelist[c-1].index].node[0] || edgelist[verticelist[t].index].node[1]==edgelist[verticelist[c-1].index].node[1]))
					{
						edgelist[   verticelist[t].index   ].color= -i;
						verticelist[t].degree--;
						w2=t;
						for (M=c; M<nedge; ++M)
							{
								if (	edgelist[verticelist[M].index].color!=-i && (M != w2) && (edgelist[verticelist[M].index].color<1) && ( edgelist[ verticelist[M].index].node[0]==edgelist[verticelist[w2].index].node[0] || edgelist[verticelist[M].index].node[0]==edgelist[verticelist[w2].index].node[1] || edgelist[verticelist[M].index].node[1]==edgelist[verticelist[w2].index].node[0] || edgelist[verticelist[M].index].node[1]==edgelist[verticelist[w2].index].node[1]))
									{
										verticelist[M].degree--;
										verticelist[M].zwischenspeicher++;
									}
							}
					}
				if (w2>0 && j!= w2 && w2>j)
					{			
						*hilfspointer =  verticelist[j];
						verticelist[j]	=	verticelist[w2];
						verticelist[w2] =	*hilfspointer;	
						j++;
					}
				else if (w2>0 && j== w2)
					{
							j++;		
					}
			}

//for (tester1= 0 ; tester1<nedge; ++tester1)    
//			{
//				printf("//ENDE// Kante %d: index %d    Grad %d    Farbe %d   TEMP %d\n", tester1, verticelist[tester1].index, verticelist[tester1].degree, edgelist[ verticelist[tester1].index].color,verticelist[tester1].zwischenspeicher);	
//			}	printf("\n\n");   

	}
	i++;
	j= c;	
	for (tester1= j ; tester1<nedge; ++tester1)    
		{
			if (edgelist[ verticelist[tester1].index].color<1)  
				{
					verticelist[tester1].degree=verticelist[tester1].degree
																				+verticelist[tester1].zwischenspeicher;
					verticelist[tester1].zwischenspeicher=0;
				}
		}
}


zeit= stopTimer ();
printf(" Zeit ist %f \n\n" , zeit);
*/



Farbanz=i-1;  printf("Farbanzahl ist %d\n", Farbanz); 
break;

case 6:  ///NTL////



printf("/////////////NTL//////////////\n\n");
printf("Anzahl der Kanten: %d\n\n", nedge);


startTimer ();

//Kantenarray initialisieren
for (j=0; j< (2*nedge) ; ++j)
  {
		Kantenarray[ j ].Kante=0;
		Kantenarray[ j ].Knoten=0;
	}


/*
//nur zum testen da
printf(" Anfangspos: ");	
for (t=0; t < nvertice  ; ++t)
	{	
		printf("%d ", Anfangspos[t]);	
	}

//nur zum testen da
printf(" aktpos: ");	
for (t=0; t < nvertice  ; ++t)
	{	
		printf("%d ", aktpos[t]);	
	}
printf("\n\n" );
//nur zum testen da
printf(" hilfsarray: ");	
for (t=0; t < nvertice  ; ++t)
	{	
		printf("%d ", hilfsarray[t]);	
	}
printf("\n\n" );
*/





for (j=0; j<nedge; ++j)
    {

			Knoten1=edgelist[j].node[0];
			Knoten2=edgelist[j].node[1];
			
			//fuege die kante j in die adjazenzliste ein
			if ( knotenliste[Knoten1].kanteknoten == NULL )  
				{
					knotenliste[Knoten1].kanteknoten=new_list(knotenliste,Knoten1,Knoten2,j);
				}
			else
				{
					knotenliste[Knoten1].kanteknoten=insert(knotenliste,Knoten1,Knoten2,j);
				}

			if ( knotenliste[Knoten2].kanteknoten == NULL )  
				{
					knotenliste[Knoten2].kanteknoten=new_list(knotenliste,Knoten2,Knoten1,j);
				}
			else
				{
					knotenliste[Knoten2].kanteknoten=insert(knotenliste,Knoten2,Knoten1,j);
				}

			//fuege die Kante in die Tabelle ein (2 mal)
			Kantenarray[ aktpos[ Knoten1 ]-1 ].Kante= j;
			Kantenarray[ aktpos[ Knoten1 ]-1 ].Knoten= Knoten2;
			aktpos[ Knoten1 ]++;

			Kantenarray[ aktpos[ Knoten2 ]-1 ].Kante= j;
			Kantenarray[ aktpos[ Knoten2 ]-1 ].Knoten= Knoten1;
			aktpos[ Knoten2 ]++;

			//bestimme w
			w=Knoten1;  

			// initialisiere W[ k ]
			for (k=0; k<colors; ++k)
				{	
					W[k]=-1;
					usedmiscol[ k ]=-1;
					Fan[k]=-2;
					Endknoten[k]=-1;
					miscolEndknoten[k]=-1;	
						
				}
			
			//bestimme W[ k ]
			struct node *run;
			struct node *run2;
			run= knotenliste[w].kanteknoten;
			miscolzaehler=0;

			while ( run  != NULL   )
				{	
					if( edgelist[ run->Kante ].color != 0  )
						{	
							W[ edgelist[ run->Kante ].color -1 ]= run->Kante ;
						}

					//bestimme die xi und deren miscol
					xi= run->Knoten;
					//printf(" xi ist %d\n", run->Knoten);
					color=0;
					zaehler=0;	
					
					while (knotenliste[xi].miscol == 0)
						{ 
							color++;
							zaehler=0;
		
							run2= knotenliste[xi].kanteknoten;
			
							while ( run2  != NULL   )
								{
									if( edgelist[ run2 -> Kante].color == color)
										{
											zaehler++;
										}
									run2=run2->right;
								}
							if(zaehler == 0)
								{
									usedmiscol[ miscolzaehler ]= xi;
									//printf("miscolzaehler ist %d\n", miscolzaehler);
									miscolzaehler++;								
									knotenliste[xi].miscol=color;
									//if(color>q)
									//		{
									//			printf("ERROR!!!! Miscolor wird zu hoch gewaehlt!\n");
									//		}
								}
						}
					//gehe zu der naechsten Kante
					run=run->right;
				}

			//bestimme miscol von w
			int miscolw;
			miscolw = 0;
			for (k=0; k<colors; ++k)
				{	
					if( W[k] == -1 && miscolw == 0   )
						{	
							miscolw= k+1;
							break;
						}
				}
			

			//Bestimme den groesstmoeglichen Fan
			s=0;
			v=Knoten2;
		
			Fan[0]=j;	
			Endknoten[0]= v;
			miscolEndknoten[0]= knotenliste[v].miscol;

			for (t=1; t < colors  ; ++t)
				{	//printf("%d und %d : ", W[knotenliste[v].miscol-1], knotenliste[v].miscol );								
					if( W[  knotenliste[v].miscol -1 ] >  -1 )
						{	
							s++;																		//kann man vllt auch weglassen
							Fan[t]= W[ knotenliste[v].miscol -1 ];

							W[ knotenliste[v].miscol -1 ] = -2;  		
												
							if( edgelist[  Fan[t] ].node[0] == w )   //geht das auch anders?
								{	
									v=edgelist[  Fan[t]   ].node[1];
									Endknoten[t]= v;
									miscolEndknoten[t]=  knotenliste[v].miscol ;
								}
							else
								{	
									v=edgelist[ Fan[t]   ].node[0];
									Endknoten[t]= v;
									miscolEndknoten[t]=  knotenliste[v].miscol ;
								}
						}
					else 
						{ 
							break;
						}
				}

			//Bestimme ws
			ws= Endknoten[ s ];

			/*
			//nur zum testen da
			printf("w ist %d, s ist %d, ws ist %d\n", w, s , ws);
			//nur zum testen da
			printf("%d.W[k]: ", j );
			for (k=0; k<colors; ++k)
			{	
				printf("%d ", W[k] );
			}
			printf("\n");
			
			//nur zum testen da
			printf("%d.Fan: ", j );
			for (k=0; k<colors; ++k)
			{	
				printf("%d ", Fan[k] );
			}
			printf("\n");
			printf(" xs ist %d\n", ws);
			
			
			//nur zum testen da
			printf("%d.Endknoten: ", j );
			for (k=0; k<colors; ++k)
			{	
				printf("%d ", Endknoten[k] );
			}
			printf("\n");
			
			//nur zum testen da
			printf("%d.miscolEndknoten: ", j );
			for (k=0; k<colors; ++k)
			{	
				printf("%d ", miscolEndknoten[k] );
			}
			printf("\n");

			//nur zum testen da
			printf("%d.miscols: ", j );
			for (t=0; t < nvertice  ; ++t)
			{	
				printf("%d ", knotenliste[t].miscol );
			}
			printf("\n\n\n");

			//nur zum testen da
			printf("%d.usedmiscol: ", j );
			for (t=0; t < colors  ; ++t)
			{	
				printf("%d ",usedmiscol[t]);
			}
			printf("\n\n\n");
			*/

		int wsmiscol =	knotenliste[ws].miscol;
				//printf(" miscol(xs) ist %d\n", wsmiscol);
		//Case 1 im Alg
		if( W[ wsmiscol -1 ] == -1 )								
			{			//printf("%d BIN IN CASE1\n",j);
				
				//geben den (w,xi) die Farbe miscol(xi)
				for (t=0; t < s+1  ; ++t)
					{
						edgelist[ Fan[t] ].color = miscolEndknoten[t];
					}
			}

		//Case 2 im Alg
		else
			{			//printf("%d BIN IN CASE2\n",j);

				//berechne xjmin1
				int xjmin1=0;
				for (r=0; r < s+1  ; ++r)
					{	
						if(  miscolEndknoten[ r ]	 == wsmiscol  )
							{
								xjmin1 = Endknoten[ r ];	
								break;
							}	
					}

				//nur zum testen da
				//printf(" xjmin1 ist %d\n",xjmin1);
			

				//benutze miscolw und miscolws fuer den Pfad
				c=0;
				int start;
				start = xjmin1;
				
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						//bestimme akt. Farbe von den alternierenden
						if( c%2==0  )
							{
								Farbe= miscolw;
							}
						else
							{
								Farbe= wsmiscol;
							}
									//printf("///Farbe ist %d\n", Farbe);
						i=0;
								//printf("///Schleife von %d bis %d\n", Anfangspos[ start ]-1, aktpos[ start ]-2);
						for (t= Anfangspos[ start ]-1 ; t < aktpos[ start ]-1  ; ++t)  //evtl anpassen
							{
			
//printf("///Bin bei Kante %d mit Farbe %d\n", Kantenarray[ t ], edgelist[  Kantenarray[ t ] ].color);
								if( edgelist[  Kantenarray[ t ].Kante ].color == Farbe  )
									{			//printf("///Bin im Prozess\n");
										c++;
										i++;
										Pfadarray[ k ]= Kantenarray[ t ].Kante;	


										start =  Kantenarray[ t ].Knoten;  

										/*
										if( edgelist[ Kantenarray[ t ].Kante ].node[0] == start  )
											{//	printf("///Bin im if\n");
												start = edgelist[ Kantenarray[ t ].Kante ].node[1] ;
											}
										else
											{	//printf("///Bin im else\n");
												start = edgelist[ Kantenarray[ t ].Kante ].node[0] ;
											}
										*/

										break;
									}	
							}
						if( i == 0  )
							{
								break;
							}
					}

				/*
				//nur zum testen da
				printf("Q: Pfadarray in Farben %d und %d von %d: ", miscolw, wsmiscol, xjmin1);
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						printf("%d ", Pfadarray[ k ]);
					}printf("c ist %d\n\n", c); */
				

				//wird w erreicht?
				int werreicht;
				werreicht=0;
				
				if( start == w  )
					{
						werreicht=1;
					}

				//falls Pfad Q nicht in w endet
				if ( werreicht == 0)	
					{	//printf("%d BIN IN CASE 2.1\n",j);	
						//printf("Pfad Q endet nicht in w!\n");
						
						//Pfad leer?
						if ( c != 0)	
							{			//printf("Q: Es muss getauscht werden\n");


						/*		
				//nur zum testen da
				printf("Q: Pfadarray in Farben %d und %d: \n", miscolw, wsmiscol);
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						printf("%d Farbe %d\n ", Pfadarray[ k ], edgelist[ Pfadarray[ k ] ].color);
					}printf("c ist %d\n\n", c);*/
				


								//Farbtausch	
								for (k= 0 ; k < nvertice  ; ++k)  
									{	
										if( Pfadarray[ k ] != -1 )
											{
													
												if( edgelist[  Pfadarray[ k ] ].color == miscolw  )
													{
														edgelist[  Pfadarray[ k ] ].color = wsmiscol;
													}
												else
													{
														edgelist[  Pfadarray[ k ] ].color = miscolw;
													}
											}
										else
											{
												break;
											}
									}

				/*
				//nur zum testen da
				printf("Q: Pfadarray in Farben %d und %d: \n", miscolw, wsmiscol);
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						printf("%d Farbe %d\n ", Pfadarray[ k ], edgelist[ Pfadarray[ k ] ].color);
					}printf("c ist %d\n\n", c);*/
				


							}
						else
							{
								//printf("Q: kein Tausch noetig!\n");
							}
	
						//Faerbe um			
						for (t=0; t < r  ; ++t)
							{
								edgelist[ Fan[t] ].color = miscolEndknoten[t];
							}
						edgelist[ Fan[r] ].color = miscolw;


				//Pfadarray neu initialisieren
						for (t=0; t < c  ; ++t)
							{
								Pfadarray[ t ] = -1;
							}		



					}

				//falls Pfad Q in w endet, endet R nicht in w!
				else
					{	//printf("%d BIN IN CASE 2.2\n",j);			
						//pfadarray neu initialisieren														
						for (t=0; t < c  ; ++t)
							{
								Pfadarray[ t ]=-1;
							}

						//benutze miscolw und miscolws fuer den Pfad
						c=0;
						int start;
						start = ws;
				
						for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
							{	
								//bestimme akt. Farbe von den alternierenden
								if( c%2==0  )
									{
										Farbe= miscolw;
									}
								else
									{
										Farbe= wsmiscol;
									}

								i=0;
								for (t= Anfangspos[ start ]-1 ; t < aktpos[ start ]-1  ; ++t)  //evtl anpassen
									{
										if( edgelist[  Kantenarray[ t ].Kante ].color == Farbe  )
											{
												c++;
												i++;
												Pfadarray[ k ]= Kantenarray[ t ].Kante;	

												start = Kantenarray[ t ].Knoten;

												/*
												if( edgelist[ Kantenarray[ t ].Kante ].node[0] == start  )
													{
														start = edgelist[ Kantenarray[ t ].Kante ].node[1];
													}
												else
													{
														start = edgelist[ Kantenarray[ t ].Kante ].node[0];
													}
												*/

												break;
											}	
									}
								if( i == 0  )
									{
										break;
									}
							}

						/*
						//nur zum testen da
						printf("R: Pfadarray in Farben %d und %d von %d: ", miscolw, wsmiscol, ws);
						for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
							{	
								printf("%d ", Pfadarray[ k ]);
							}printf("c ist %d\n\n", c);
						*/


						//wird w erreicht?
						int werreicht;
						werreicht=0;
				
						if( start == w  )
							{
								werreicht=1;
							}

						//falls Pfad R nicht in w endet
						if ( werreicht == 0)	
							{		//printf("Pfad R endet nicht in w!\n");
									
								//Pfad leer?
								if ( c != 0)	
									{	//printf("R: Es muss getauscht werden\n");
					
				/*
				//nur zum testen da
				printf("R: Pfadarray in Farben %d und %d: \n", miscolw, wsmiscol);
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						printf("%d Farbe %d\n ", Pfadarray[ k ], edgelist[ Pfadarray[ k ] ].color);
					}printf("c ist %d\n\n", c);*/
										
										//Farbtausch	
										for (k= 0 ; k < nvertice  ; ++k)  
											{	
												if( Pfadarray[ k ] != -1 )
													{
													
														if( edgelist[  Pfadarray[ k ] ].color == miscolw  )
															{	
																edgelist[  Pfadarray[ k ] ].color = wsmiscol;
															}
														else
															{
																edgelist[  Pfadarray[ k ] ].color = miscolw;
															}
													}
												else
													{
														break;
													}
											}


				/*
				//nur zum testen da
				printf("R: Pfadarray in Farben %d und %d: \n", miscolw, wsmiscol);
				for (k= 0 ; k < nvertice  ; ++k)  //while schleife?
					{	
						printf("%d Farbe %d\n ", Pfadarray[ k ], edgelist[ Pfadarray[ k ] ].color);
					}printf("c ist %d\n\n", c);
					*/

									}
								else
									{
										//printf("R: kein Tausch noetig!\n");
									}
	
								//Faerbe um
								for (t=0; t < s  ; ++t)
									{
										edgelist[ Fan[t] ].color = miscolEndknoten[t];
									}
								edgelist[ Fan[s] ].color = miscolw;

						//pfadarray neu initialisieren														
						for (t=0; t < c  ; ++t)
							{
								Pfadarray[ t ]=-1;
							}



							}
					}
			}// ende von case 2

				
		// W[ ] und Fan[ ] neu initialsieren
		for (tester1=0; tester1<q; ++tester1)
			{	
				W[tester1]=-1;
				Fan[tester1]=-2;
			}


		//miscol neu initialisieren									
		for (t=0; t < colors  ; ++t)
			{	
				if( usedmiscol[ t ] != -1  )
					{
						knotenliste[ usedmiscol[ t ] ].miscol = 0;
						usedmiscol[ t ]= -1;
					}
				else
					{
						break;
					}
			}	 

		/*
		//nur zum testen da
			printf("%d.miscols: ", j );
			for (t=0; t < nvertice  ; ++t)
			{	
				printf("%d ", knotenliste[t].miscol );
			}
			printf("\n\n\n");

		//nur zum testen da
			printf("%d.usedmiscol: ", j );
			for (t=0; t < colors  ; ++t)
			{	
				printf("%d ",usedmiscol[t]);
			}
			printf("\n\n\n");
		*/


		/*
		//nur zum testen da
		for (r=0; r<j+1; ++r)
   	 {
			printf("Kante %d hat Knoten %d %d und Farbe %d\n", r, edgelist[r].node[0], edgelist[r].node[1], edgelist[r].color);
 			}
		printf("\n\n");*/
		
}

/*
//nur zum testen da
printf(" Kantenarray: " );
for (j=0; j< (2*nedge) ; ++j)
  {
		printf("%d " , Kantenarray[ j ]);
	}printf("\n\n");

//nur zum testen da
printf(" aktpos: ");	
for (t=0; t < nvertice  ; ++t)
	{	
		printf("%d ", aktpos[t]);	
	}
*/

zeit= stopTimer ();
printf("\n\nBenoetigte Zeit: %f \n" , zeit);


printf("//////////////////////////\n\n");
break;

case 7: //nichts//
break;	
}
	
		
switch(postprocessing)
{
case 1: //Kanten//

for (j=0; j<nedge; ++j)
    {
		printf("Kante %d hat Knoten %d %d und Farbe %d\n", j, edgelist[j].node[0], edgelist[j].node[1], edgelist[j].color);
 		}

break;

case 2: //Knoten//

printf("\n");

	for (j=0; j<nedge; ++j)
    {
		printf("Knoten %d hat Grad %d und alten Index %d\n", j, verticelist[j].degree, verticelist[j].index);
		}

printf("\n");

for (j=0; j<nedge; ++j)
    {
			printf("Knoten %d-Nachbarn mit altem Index %d:      ", j,verticelist[j].index);
			for (nachbarzaehler=0; nachbarzaehler<verticelist[j].degree; ++nachbarzaehler)
			{
				printf(" %d ", verticelist[j].nachbar[nachbarzaehler]);
			}
		printf("\n");
		}

break;

case 3: //Anz. der Elemente in einer Farbe und Zulaessigkeit

gesamt=0;
printf("\n");
if (Algorithmus != 6)
	{
		for (j=0; j<max; ++j)
    	{
				printf("Anzahl der Elemente von Farbe %d: %d\n", j+1, Farbanzahl[j]);
				gesamt=gesamt+Farbanzahl[j];
 			}
		if( gesamt== nedge)
			{
				printf("\nAnzahl stimmt :) \n");
			}
		else printf("\nAnzahl stimmt nicht!!! :( \n");
	} 
//else
//	{	
//		int f=0;
//		printf("\n");
//		struct node *aktuelll;
//		for(k=0; k<q; ++k)  
//			{
//				aktuelll= farbenliste[k].Elemente;
//				//printf("Farbliste%d: " ,k+1);
//				if (farbenliste[k].Elemente != NULL)
//				{				
//					while (aktuelll->right != NULL)
//					{	f++;
//						//printf(" %d " ,aktuelll->data);
//						aktuelll=aktuelll->right;
//					}
//					//printf(" %d ",aktuelll->data);
//					f++;
//					printf("Anzahl der Elemente von Farbe %d: %d\n" ,k+1, f);
//					gesamt=gesamt+f;
//					f=0;
//				}		
//				else if(farbenliste[k].Elemente == NULL)
//					{	//printf(" -----------");
//						printf("Anzahl der Elemente von Farbe %d: %d\n" ,k+1, f);
//					}
//			}
//		if( gesamt== nedge)
//			{
//				printf("\nAnzahl stimmt :) \n");
//			}
//		else printf("\nAnzahl stimmt nicht!!! :( \n");
//	}
//		printf(" Anzahl der Kanten ist %d \n", nedge);
//		printf("\n");  



break;

case 4: 

//Faerbung zulaessig?

m=0;
for (j=0; j<nedge; ++j)
	{
		for (k=0; k<prueferlist[j].pruefdegree; ++k)
			{ 
				if (edgelist[j].color==edgelist[prueferlist[j].pruefnachbarn[k]].color)
					{
						printf(" Farbe von Kante %d und Kante %d ist gleich!\t Naemlich: %d! \n", j, prueferlist[j].pruefnachbarn[k],edgelist[j].color);
						m++;
					}			
			}//printf("\n");
	}

if (m==0)
	{
		printf("\n\nFaerbung zulaessig =) \n");
		m++;
	}
else
	{
		printf("\n\nFaerbung ist unzulaessig!!!\n");
	}

break;

case 5:

printf("\n");
struct node *look;
for (j=0; j<nvertice; ++j)
	{	
		look= knotenliste[j].kanteknoten;
		printf("Knoten %d: ",j );
		while ( look!= NULL   )
		{	printf("Kante %d - %d // " ,look->Kante, look->Knoten);	
			look=look->right;
		}printf("\n");
	}

break;
} 

/*
//Speicher freen
struct node *temp2;

if( (Algorithmus !=4 && Algorithmus !=6) ||postprocessing == 4 || postprocessing == 2)													
	{
		for (j=0; j<nedge; ++j)														
			{
				free(verticelist[j].nachbar);
				free(prueferlist[j].pruefnachbarn);
			}
	}
//free(knotenliste);
free(verticelist);
free(prueferlist);
free(hilfspointer); 
*/
return Farbanz;  
}

/*
 * Wrapper routine that can be called from Fortran
 *
 * neq  :        number of equations (=vertices)
 * nedge:        number of edges
 * ncolor:       maximum number of colors on input,
 *               actual number of colors on output
 * IedgeList:    pointer to the edge List
 * IedgeListIdx: index pointer separating groups of edges with the same color
 * 
 * Example:
 *   IedgeListIdx = [1,5,10,18];
 *   IedgeList    = [(1,2),(3,4),(5,7),(8,3),...]
 *
 * Edges  1..4  have color C1,
 * edges  5..9  have color C2,
 * edges 10..17 have color C3,...
 */
void FNAME(regroupedgelist)(int *neq, int *nedge, int *ncolor,
			    int **IedgeList, int **IedgeListIdx)
{
  // Allocate list of edges in C-format
  Edge *edgelist = (Edge*) malloc(*nedge*sizeof(Edge));
  
  int iedge,icolor;
  int *d_IedgeList = (int*) IedgeList;
  int *d_IedgeListIdx = (int*) IedgeListIdx;
  
  // Fill list of edges
  for (iedge=0; iedge<(*nedge); ++iedge) {
    // Empty color
    edgelist[iedge].color = 0;
    
    // 0-based vertex numbers
    edgelist[iedge].node[0] = d_IedgeList[6*iedge]-1;
    edgelist[iedge].node[1] = d_IedgeList[6*iedge+1]-1;

    // auxiliary data (do not care about)
    edgelist[iedge].data[0] = d_IedgeList[6*iedge+2];
    edgelist[iedge].data[1] = d_IedgeList[6*iedge+3];
    edgelist[iedge].data[2] = d_IedgeList[6*iedge+4];
    edgelist[iedge].data[3] = d_IedgeList[6*iedge+5];
  }
  
  // Apply edge coloring algorithm
  int ncolors = edgeColoring(*neq, *nedge, edgelist);
  if (ncolors > *ncolor) {
    printf("Error: number of colors exceds maximum number of colors!\n");
    return;
  }
  
  // Clear index array
  for (icolor=0; icolor<(*ncolor); ++icolor)
    d_IedgeListIdx[icolor] = 0;

  // Loop over all color groups and over all edges
  int icount=0;
  for (icolor=0; icolor<ncolors; icolor++) {
    d_IedgeListIdx[icolor] = icount+1; // we are 1-based in Fortran

    for (iedge=0; iedge<(*nedge); ++iedge) {
      if (edgelist[iedge].color == icolor) {
	
	// 1-based vertex numbers
	d_IedgeList[6*icount  ] = edgelist[iedge].node[0]+1;
	d_IedgeList[6*icount+1] = edgelist[iedge].node[1]+1;

	// auxiliary data
	d_IedgeList[6*icount+2] = edgelist[iedge].data[0];
	d_IedgeList[6*icount+3] = edgelist[iedge].data[1];
	d_IedgeList[6*icount+4] = edgelist[iedge].data[2];
	d_IedgeList[6*icount+5] = edgelist[iedge].data[3];

	icount++;
      }
    }
  }
  
  // Fill unused groups
  for (icolor=ncolors; icolor<(*ncolor); ++icolor)
    d_IedgeListIdx[icolor] = icount+1;
  
  // Deallocate list of edges
  free(edgelist);
}
