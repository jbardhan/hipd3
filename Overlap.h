#ifndef __OVERLAP_H__
#define __OVERLAP_H__
#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
void Preconditioner_do_leaf_leaf_translation(SMatrix mat, Cube cube, Cube srccube);
void Solv_ecf_qual_writematlabfile(char *filename, Tree tree);
void Preconditioner_fill_overlap_get_adjacent_cubes(Cube cube, unsigned int *numadj,
                                                    Cube **adjcubes);
void Preconditioner_fill_overlap_recurse(Preconditioner preconditioner, Cube cube, Panel *panels,
                                         real idiel, real odiel);
void Preconditioner_fill_overlap_ecf_qual_cav(Preconditioner preconditioner, Tree tree,
                                              Panel* panels, unsigned int numpoints, unsigned int numpanels,
                                              real idiel, real odiel);
void Preconditioner_multiply(Vector Px, Preconditioner preconditioner, Vector x, unsigned int numpanels);
#endif
