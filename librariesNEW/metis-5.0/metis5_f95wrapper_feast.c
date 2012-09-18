/* The Intel Fortran Runtime Library does not accept a NULL pointer being
 * passed from Fortran to C when bounds checking is enabled. In such cases it
 * aborts with
 *     forrtl: severe (408): fort: (7): Attempt to use pointer ... when it
 *             is not associated with a target
 * So, calling METIS libraries with NULL arguments (to tell METIS 5.0 to use
 * defaults) is not feasible directly from Fortran 95 code.
 * As a workaround, provide here wrapper functions to METIS 5.0 with an
 * interface that omits the NULL pointers and hardwire the NULL pointers in
 * the C wrapper.
 */
#include <stdio.h>
#include "include/metis.h"

int feastmetis_partgraphrecursive(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
          idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
          idx_t *nparts, idx_t *options, idx_t *objval, idx_t *part)
{
  return METIS_PartGraphRecursive(nvtxs, ncon, xadj,
				  adjncy, vwgt, vsize, adjwgt,
				  nparts, NULL, NULL, options, objval, part);
}

int feastmetis_partgraphrecursive_(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
          idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
          idx_t *nparts, idx_t *options, idx_t *objval, idx_t *part)
{
  return feastmetis_partgraphrecursive(nvtxs, ncon, xadj,
				       adjncy, vwgt, vsize, adjwgt,
				       nparts, options, objval, part);
}

int feastmetis_partgraphrecursive__(idx_t *nvtxs, idx_t *ncon, idx_t *xadj,
          idx_t *adjncy, idx_t *vwgt, idx_t *vsize, idx_t *adjwgt,
          idx_t *nparts, idx_t *options, idx_t *objval, idx_t *part)
{
  return feastmetis_partgraphrecursive(nvtxs, ncon, xadj,
				       adjncy, vwgt, vsize, adjwgt,
				       nparts, options, objval, part);
}

int feastmetis_partgraphkway(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy,
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts,
          idx_t *options, idx_t *objval,
          idx_t *part)
{
  return METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy,
			     vwgt, vsize, adjwgt, nparts,
			     NULL, NULL, options, objval, part);
}

int feastmetis_partgraphkway_(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy,
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts,
          idx_t *options, idx_t *objval,
          idx_t *part)
{
  return feastmetis_partgraphkway(nvtxs, ncon, xadj, adjncy,
				  vwgt, vsize, adjwgt, nparts,
				  options, objval, part);
}

int feastmetis_partgraphkway__(idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy,
          idx_t *vwgt, idx_t *vsize, idx_t *adjwgt, idx_t *nparts,
          idx_t *options, idx_t *objval,
          idx_t *part)
{
  return feastmetis_partgraphkway(nvtxs, ncon, xadj, adjncy,
				  vwgt, vsize, adjwgt, nparts,
				  options, objval, part);
}
