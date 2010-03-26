

#include   <graham.h>
#include   <stdio.h>
#include   <math.h>
#include   <stdlib.h>

int n = 0;                         /* Actual # of points */
int ndel = 0;                   /* Number deld */

int     AreaSign5( tPointi a, tPointi b, tPointi c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}

int   Compare( const void *tpi, const void *tpj )
{
   int a;             /* area */
   int x, y;          /* projections of ri & rj in 1st quadrant */
   tPoint pi, pj;
   pi = (tPoint)tpi;
   pj = (tPoint)tpj;

   a = AreaSign5( P[0].v, pi->v, pj->v );
   if (a > 0)
      return -1;
   else if (a < 0)
      return 1;
   else { /* Collinear with P[0] */
      x =  abs( pi->v[X] -  P[0].v[X] ) - abs( pj->v[X] -  P[0].v[X] );
      y =  abs( pi->v[Y] -  P[0].v[Y] ) - abs( pj->v[Y] -  P[0].v[Y] );

      ndel++;
      if ( (x < 0) || (y < 0) ) {
         pi->del = TRUE;
         return -1;
      }
      else if ( (x > 0) || (y > 0) ) {
         pj->del = TRUE;
         return 1;
      }
      else { /* points are coincident */
         if (pi->vnum > pj->vnum)
             pj->del = TRUE;
         else 
             pi->del = TRUE;
         return 0;
      }
   }
}



std::vector<int> graham::compConvexHull()
{

   std::vector<int> indices;

   tStack   top;

   n = ReadPoints();
   FindLowest();
   PrintPoints();
   qsort(
      &P[1],             /* pointer to 1st elem */
      n-1,               /* number of elems */
      sizeof( tsPoint ), /* size of each elem */
      Compare            /* -1,0,+1 compare function */
   );
  // printf("After sorting, ndel = %d:\n", ndel);
   //PrintPoints();
   if (ndel > 0) {
      Squash();
      printf("After squashing:\n");
   }

   top = Graham();
   //printf("Hull:\n");
   
   indices = PrintStack( top );
   
   return indices;
   
}

graham::graham()
{
}

void   graham::getPoints(VECTOR2* pvPts, int num, int scale)
{
	m_inumPoints = num;
	sc = scale;
	m_vPoints = new VECTOR2[num];
	for(int i = 0; i < num; i++)
	{
		m_vPoints[i] = pvPts[i];
	}
	
}

void    graham::getPoints(vector<VECTOR2> vPts, int scale)
{
	m_inumPoints = (int)vPts.size();
	sc = scale;
	m_vPoints = new VECTOR2[m_inumPoints];
	for(int i = 0; i < m_inumPoints; i++)
	{
		m_vPoints[i] = vPts[i];
	}
	
}

graham::~graham()
{
	if(m_vPoints)
	{
		delete[] m_vPoints;
		m_vPoints = NULL;
	}
}

void   graham::FindLowest( void )
{
   int i;
   int m = 0;   /* Index of lowest so far. */

   for ( i = 1; i < n; i++ )
      if ( (P[i].v[Y] <  P[m].v[Y]) ||
          ((P[i].v[Y] == P[m].v[Y]) && (P[i].v[X] > P[m].v[X])) ) 
         m = i;
  
   Swap(0,m); /* Swap P[0] and P[m] */
}

void	graham::Swap( int i, int j )
{
   int temp;
   bool tmp;
   /* Uses swap macro. */

   SWAP( temp, P[i].vnum,   P[j].vnum );
   SWAP( temp, P[i].v[X],   P[j].v[X] );
   SWAP( temp, P[i].v[Y],   P[j].v[Y] );
   SWAP( tmp, P[i].del, P[j].del );

}
/*---------------------------------------------------------------------
Compare: returns -1,0,+1 if p1 < p2, =, or > respectively;
here "<" means smaller angle.  Follows the conventions of qsort.
---------------------------------------------------------------------*/



/*---------------------------------------------------------------------
Pops off top elment of stack s, frees up the cell, and returns new top.
---------------------------------------------------------------------*/
tStack   graham::Pop( tStack s )
{
   tStack top;

   top = s->next;
   FREE( s );
   return top;
}

/*---------------------------------------------------------------------
Get a new cell, fill it with p, and push it onto the stack.
Return pointer to new stack top.
---------------------------------------------------------------------*/
tStack   graham::Push( tPoint p, tStack top )
{
   tStack   s;

   /* Get new cell and fill it with point. */
   NEW( s, tsStack );
   s->p = p;
   s->next = top;
   return s;
}
/*---------------------------------------------------------------------
---------------------------------------------------------------------*/
std::vector<int>   graham::PrintStack( tStack t )
{
	std::vector<int> ind;

   if (!t) printf("Empty stack\n");
   while (t) { 
      //printf("vnum=%d\tx=%d\ty=%d\n", 
      //       t->p->vnum,t->p->v[X],t->p->v[Y]);
	  ind.push_back(t->p->vnum);
      tStack top = t->next;
      FREE(t);
      t = top;
	  
   }
   return ind;
}

/*---------------------------------------------------------------------
Performs the Graham scan on an array of angularly sorted points P.
---------------------------------------------------------------------*/
tStack   graham::Graham()
{
   tStack   top;
   int i;
   tPoint p1, p2;  /* Top two points on stack. */

   /* Initialize stack. */
   top = NULL;
   top = Push ( &P[0], top );
   top = Push ( &P[1], top );

   /* Bottom two elements will never be removed. */
   i = 2;

   while ( i < n ) {
      //printf("Stack at top of while loop, i=%d, vnum=%d:\n", i, P[i].vnum);
      //PrintStack( top );
      if( !top->next) printf("Error\n"),exit(EXIT_FAILURE);
      p1 = top->next->p;
      p2 = top->p;
      if ( Left( p1->v , p2->v, P[i].v ) ) {
         top = Push ( &P[i], top );
         i++;
      } else    
         top = Pop( top );
      //printf("Stack at bot of while loop, i=%d, vnum=%d:\n", i, P[i].vnum);
      //PrintStack( top );
      //putchar('\n');
   }

   return top;

}

/*---------------------------------------------------------------------
Squash removes all elements from P marked del.
---------------------------------------------------------------------*/
void   graham::Squash( void )
{
   int i, j;
   i = 0; j = 0;
   /*printf("Squash: n=%d\n",n);*/
   while ( i < n ) {
      /*printf("i=%d,j=%d\n",i,j);*/
      if ( !P[i].del ) { /* if not marked for deletion */
         Copy( i, j ); /* Copy P[i] to P[j]. */
         j++;
      }
      /* else do nothing: del by skipping. */
      i++;
   }
   n = j;
   
   PrintPoints();
}

void	graham::Copy( int i, int j )
{
   P[j].v[X] = P[i].v[X];
   P[j].v[Y] = P[i].v[Y];
   P[j].vnum = P[i].vnum;
   P[j].del = P[i].del;
}

int* graham::cpyPoints(int num, tStack t)
{

	int* indices = new int[num];
	int j = 0;
	if (!t) printf("Empty stack\n");
	
	while (t)
	{ 
			 printf("vnum=%d\tx=%d\ty=%d\n", 
			 t->p->vnum,t->p->v[X],t->p->v[Y]); 
			 indices [j] = t->p->vnum  ;
			 t = t->next;
			 j++;
		
	}

	return indices;


}

/*---------------------------------------------------------------------
Returns twice the signed area of the triangle determined by a,b,c.
The area is positive if a,b,c are oriented ccw, negative if cw,
and zero if the points are collinear.
---------------------------------------------------------------------*/
int     graham::Area2( tPointi a, tPointi b, tPointi c )
{
   return
          (b[X] - a[X]) * (c[Y] - a[Y]) -
          (c[X] - a[X]) * (b[Y] - a[Y]);
}

/*---------------------------------------------------------------------
Returns true iff c is strictly to the left of the directed
line through a to b.
---------------------------------------------------------------------*/
bool   graham::Left( tPointi a, tPointi b, tPointi c )
{ 
   return  Area2( a, b, c ) > 0;
}

/*---------------------------------------------------------------------
Reads in the coordinates of the points from stdin,
puts them into P, and returns n, the number of vertices.
Initializes other fields of point structure.
---------------------------------------------------------------------*/
int    graham::ReadPoints(  )
{

	int n = 0;
	int i;

   for(i = 0; i < m_inumPoints; i++)
   {
	   
       double x = m_vPoints[i].x * sc;
       double y = m_vPoints[i].y * sc;
	   /*P[i].v[0] =(int)m_vPoints[i].x;
	   P[i].v[1] =(int)m_vPoints[i].y;*/
	   P[i].v[0] =(int)x;
	   P[i].v[1] =(int)y;
	   P[i].vnum = i;
	   P[i].del = false;
   }

   n = i;

  // while ( (n < PMAX) && (scanf("%d %d",&P[n].v[0],&P[n].v[1]) != EOF) ) {
  //    P[n].vnum = n;
  //    P[n].del = FALSE;
  //    /*printf("vnum=%3d, x=%4d, y=%4d, del=%d\n", 
	 //P[n].vnum, P[n].v[X], P[n].v[Y], P[n].del);*/
  //    ++n;
  // }
   if(n >= PMAX)
   {
     printf("Error in ReadPoints:  too many points; max is %d\n", PMAX);
      exit(EXIT_FAILURE);
   }//end if

   return n;
}

void   graham::PrintPoints( void )
{
   int   i;

   //printf("Points:\n");
   for( i = 0; i < n; i++ );
      //printf("vnum=%3d, x=%4d, y=%4d, del=%d\n", 
	    // P[i].vnum, P[i].v[X], P[i].v[Y], P[i].del);
}

void   graham::PrintPostscript( tStack t)
{
   int   i;
   int xmin, ymin, xmax, ymax;

   xmin = xmax = P[0].v[X];
   ymin = ymax = P[0].v[Y];
   for (i = 1; i < n; i++) {
      if      ( P[i].v[X] > xmax ) xmax = P[i].v[X];
      else if ( P[i].v[X] < xmin ) xmin = P[i].v[X];
      if      ( P[i].v[Y] > ymax ) ymax = P[i].v[Y];
      else if ( P[i].v[Y] < ymin ) ymin = P[i].v[Y];
   }

   /* PostScript header */
   printf("%%!PS\n");
   printf("%%%%Creator: graham.c (Joseph O'Rourke)\n");
   printf("%%%%BoundingBox: %d %d %d %d\n", xmin, ymin, xmax, ymax);
   printf("%%%%EndComments\n");
   printf(".00 .00 setlinewidth\n");
   printf("%d %d translate\n", -xmin+72, -ymin+72 );
   /* The +72 shifts the figure one inch from the lower left corner */

   /* Draw the points as little circles. */
   printf("newpath\n");
   printf("\n%%Points:\n");
   for (i = 0; i < n; i++)
      printf("%d\t%d\t1  0  360\tarc\tstroke\n", P[i].v[X], P[i].v[Y]);
   printf("closepath\n");

   /* Draw the polygon. */
   printf("\n%%Hull:\n");
   printf("newpath\n");
   printf("%d\t%d\tmoveto\n", t->p->v[X], t->p->v[Y]);
   while (t) {
      printf("%d\t%d\tlineto\n", t->p->v[X], t->p->v[Y]);
      t = t->next;
   }
   printf("closepath stroke\n");
   printf("showpage\n%%%%EOF\n");
}
int     graham::AreaSign( tPointi a, tPointi b, tPointi c )
{
    double area2;

    area2 = ( b[0] - a[0] ) * (double)( c[1] - a[1] ) -
            ( c[0] - a[0] ) * (double)( b[1] - a[1] );

    /* The area should be an integer. */
    if      ( area2 >  0.5 ) return  1;
    else if ( area2 < -0.5 ) return -1;
    else                     return  0;
}
