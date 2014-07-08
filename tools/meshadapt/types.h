/*
 * This header file translates all derived types from the Featflow2
 * kernel as C structs.
 *
 */

// fparser.f90
struct t_fparserComponent {
  int ibytecodeSize;
  int iimmedSize;
  int istackSize;
  int istackPtr;
  bool buseDegreeConversion;
  bool bisVectorizable;
  short*  IbyteCode;
  double* Dimmed;
};

struct t_fparser {
  struct t_fparserComponent* Rcomp;
  char* ScompName;
  int ncomp;
  int nncomp;
};

// boundary.f90
struct t_boundary {
  int iboundarycount_g;
  int iboundarycount;
  int h_DmaxPar;
  int h_IsegCount;
  int h_Idbldatavec_handles;
  int h_Iintdatavec_handles;
  int h_Iintdatavec_fparser;
  struct t_fparser* p_rfparser;
};

// quadtreeDP.f90
struct t_quadtreeDP {
  int NFREE;
  int NDATA;
  int NVT;
  int NNVT;
  int NNODE;
  int NNNODE;
  int NRESIZE;
  double dfactor;
  int h_Data;
  int h_BdBox;
  double* p_Data;
  double* p_BdBox;
  int h_Knode;
  int* p_Knode;
};

// octreeDP.f90
struct t_octreeDP {
  int NFREE;
  int NDATA;
  int NVT;
  int NNVT;
  int NNODE;
  int NNNODE;
  int NRESIZE;
  double dfactor;
  int h_Data;
  int h_BdBox;
  double* p_Data;
  double* p_BdBox;
  int h_Knode;
  int* p_Knode;
};

// mapInt_DP.f90
struct t_mapInt_DP {
  int NA;
  int NNA;
  int NNA0;
  int nresize;
  double dfactor;
  int h_Key;
  int* p_Key;
  int* p_Kbal;
  int* p_Kparent;
  int* p_Kchild;
  int isizeData;
  int h_Data;
  double* p_Data;
};

// arraylistInt.f90
struct t_arraylistInt {
  int ntable;
  int nntable;
  int nntable0;
  int NA;
  int NNA;
  int NNA0;
  int nresize;
  double dfactor;
  int h_Ktable;
  int* p_Ktable;
  int h_Knext;
  int h_Kprev;
  int* p_Knext;
  int* p_Kprev;
  int h_Key;
  int* p_Key;
};

// hadaptaux.f90
struct t_hadapt {
  int iSpec;
  int iduplicationFlag;
  int iadaptationStrategy;
  int nsubdividemax;
  int nRefinementSteps;
  int nCoarseningSteps;
  int nSmoothingSteps;
  int drefinementTolerance;
  int dcoarseningTolerance;
  int ndim;
  int NVT0;
  int NVT;
  int increaseNVT;
  int NVBD0;
  int NVBD;
  int NBCT;
  int NEL0;
  int NEL;
  int NELMAX;
  int nGreenElements;
  int InelOfType[8];
  int InelOfType0[8];
  int h_Imarker;
  int h_IvertexAge;
  int* p_IvertexAge;
  int h_InodalProperty;
  int* p_InodalProperty;
  int h_IverticesAtElement;
  int* p_IverticesAtElement;
  int h_IneighboursAtElement;
  int* p_IneighboursAtElement;
  int h_ImidneighboursAtElement;
  int* p_ImidneighboursAtElement;
  int h_DvertexCoords1D;
  struct t_quadtreeDP rVertexCoordinates2D;
  struct t_octreeDP rVertexCoordinates3D;
  struct t_mapInt_DP* rBoundary;
  struct t_arraylistInt rElementsAtVertex;
};

// triangulation.f90
struct t_triangulation {
  int iduplicationFlag;
  int ndim;
  int NVT;
  int NMT;
  int NAT;
  int NEL;
  int NBCT;
  int NblindBCT;
  int NVBD;
  int NABD;
  int NMBD;
  int NNVE;
  int NNEE;
  int NNAE;
  int NNVA;
  int NNelAtVertex;
  int NNelAtEdge;
  int InelOfType[8];
  int nverticesPerEdge;
  int nVerticesOnAllEdges;
  int nverticesInEachElement;
  int nverticesInAllElements;
  int nadditionalVertices;
  double DboundingBoxMin[3];
  double DboundingBoxMax[3];
  int h_DvertexCoords;
  int h_IverticesAtElement;
  int h_IedgesAtElement;
  int h_IneighboursAtElement;
  int h_IelementsAtEdge;
  int h_IverticesAtEdge;
  int h_InodalProperty;
  int h_ImacroNodalProperty;
  int h_DelementVolume;
  int h_IelementsAtVertexIdx;
  int h_IelementsAtVertex;
  int h_IrefinementPatchIdx;
  int h_IrefinementPatch;
  int h_IcoarseGridElement;
  int h_IboundaryCpIdx;
  int h_IverticesAtBoundary;
  int h_IboundaryCpEdgesIdx;
  int h_IedgesAtBoundary;
  int h_IboundaryCpFacesIdx;
  int h_IfacesAtBoundary;
  int h_IelementsAtBoundary;
  int h_DvertexParameterValue;
  int h_DedgeParameterValue;
  int h_IboundaryVertexPos;
  int h_IboundaryEdgePos;
  int h_DfreeVertexCoordinates;
  int h_IelementsAtEdgeIdx3D;
  int h_IelementsAtEdge3D;
  int h_IfacesAtEdgeIdx;
  int h_IfacesAtEdge;
  int h_IfacesAtVertexIdx;
  int h_IfacesAtVertex;
  int h_IfacesAtElement;
  int h_IverticesAtFace;
  int h_IelementsAtFace;
  int h_IedgesAtFace;
  int h_IedgesAtVertexIdx;
  int h_IedgesAtVertex;
  int h_ItwistIndex;
};
