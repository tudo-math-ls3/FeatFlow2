#define EXIT_FAILURE 1
char *malloc();

#define NEW(p, type) if ((p=(type *) malloc (sizeof(type))) == NULL) {\
      printf ("NEW: Out of Memory!\n");\
      exit(EXIT_FAILURE);\
    }

#define FREE(p)  if (p) { free ((char *) p); p = NULL; }

#define SWAP(t,x,y)     { t = x; x = y; y = t; }
