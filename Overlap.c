#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Overlap.h"

void Solv_ecf_qual_writematlabfile(char *filename, Tree tree) {
   unsigned int i,j;
   Vector x = Vector_allocate(tree->numpanels);
   Vector ans = Vector_allocate(tree->numpoints);
   FILE* file = NULL;

   file = fopen(filename, "w");

   fprintf(file, "A = zeros(%u, %u);\n", tree->numpanels, tree->numpoints);

   fprintf(file, "A = [\n");

   for (i = 0; i < tree->numpanels; i++) {
      Vector_zero(x, tree->numpanels);
      x[i] = 1.0;
      Solv_ecf_multiply_qual(ans, tree, x, 4, 80);
      for (j = 0; j < tree->numpoints; j++)
         fprintf(file, "%e ", ans[j]);
      fprintf(file, "\n");
   }

   fprintf(file, "]';\n");

   Vector_free(x);
   Vector_free(ans);
}

void Preconditioner_fill_overlap_get_adjacent_cubes(Cube cube, unsigned int *numadj,
                                                    Cube **adjcubes) {
   unsigned int i, j, delta;
   unsigned int dx, dy, dz;
   
   *numadj = 1;
   for (i = 0; i < cube->numlocalcubes; i++) {
      dx = abs(cube->indices[0] - cube->localcubes[i]->indices[0]);
      dy = abs(cube->indices[1] - cube->localcubes[i]->indices[1]);
      dz = abs(cube->indices[2] - cube->localcubes[i]->indices[2]);
      if (( dx <= 1 ) && ( dy <= 1) && (dz <= 1))
         *numadj = *numadj + 1;
   }
   *adjcubes = (Cube *)calloc(*numadj, sizeof(Cube));

   (*adjcubes)[0] = cube;
   j = 1;
   for (i = 0; i < cube->numlocalcubes; i++) {
      dx = abs(cube->indices[0] - cube->localcubes[i]->indices[0]);
      dy = abs(cube->indices[1] - cube->localcubes[i]->indices[1]);
      dz = abs(cube->indices[2] - cube->localcubes[i]->indices[2]);
      if (( dx <= 1 ) && ( dy <= 1) && (dz <= 1)) {
         (*adjcubes)[j] = cube->localcubes[i];
         j++;
      }
   }
}

void Preconditioner_fill_overlap_recurse(Preconditioner preconditioner, Cube cube, Panel *panels,
                                         real idiel, real odiel) {
   unsigned int cx, cy, cz;
   SMatrix Psub, Psubinv;
   unsigned int numlocalpoints = 0;
   unsigned int numlocalpanels = 0;
   unsigned int numalllocal;
   Cube *alllocalcubes;
   unsigned int i, j, curpanel;
   unsigned int curStartCol, curStartRow, curcol;
   unsigned int localSrc, localDest;
   sreal *diag, *areas;
   Cube srccube, destcube;
   
   if (cube->leaf) {
      Preconditioner_fill_overlap_get_adjacent_cubes(cube, &numalllocal, &(alllocalcubes));
      for (i = 0; i < numalllocal; i++) {
         numlocalpoints += alllocalcubes[i]->numpointindices;
         numlocalpanels += alllocalcubes[i]->numpanelindices;
      }
      Psub    = SMatrix_allocate(numlocalpoints, numlocalpanels);
      Psubinv = SMatrix_allocate(numlocalpoints, numlocalpanels);
      diag = (sreal *)calloc(numlocalpoints, sizeof(sreal));
      areas = (sreal *)calloc(numlocalpanels, sizeof(sreal));
      curpanel = 0;
      for (i = 0; i < numalllocal; i++) {
         for (j = 0; j < alllocalcubes[i]->numpanelindices; j++)
            areas[curpanel++] = (sreal) panels[alllocalcubes[i]->panelindices[j]]->area;
      }
      // fill submatrix Psub, of P for all localcubes interacting with all local cubes
      curStartRow = 0;
		//		printf("numalllocal = %d\n", numalllocal);
      for (localDest = 0; localDest < numalllocal; localDest++) {
         destcube = alllocalcubes[localDest];
         curStartCol = 0;
         for (localSrc = 0; localSrc < numalllocal; localSrc++) {
            srccube = alllocalcubes[localSrc];
            curcol = destcube->numpanelindices;
            for (i = 0; i < destcube->numlocalcubes; i++) {
               if (srccube == destcube->localcubes[i]) {
                  SMatrix_copypiece(Psub, curStartRow, curStartCol,
                                   destcube->D_double, 0, curcol,
                                   destcube->numpointindices, srccube->numpanelindices);
	               }
               curcol += destcube->localcubes[i]->numpanelindices;
            }

				// if we've gotten this far, destcube and srccube are "well separated" via local = 1;
				for (i = 0; i < destcube->numinteractingcubes; i++) {
				  if (srccube == destcube->interactingcubes[i]) {
					 SMatrix nearbyCompressedMat = SMatrix_allocate(destcube->numpointindices, srccube->numpanelindices);
					 Preconditioner_do_leaf_leaf_translation(nearbyCompressedMat, destcube, srccube);
					 SMatrix_copypiece(Psub, curStartRow, curStartCol,
											 nearbyCompressedMat, 0, 0,
											 destcube->numpointindices, srccube->numpanelindices);
					 SMatrix_free(nearbyCompressedMat);
				  }
				}

            curStartCol += srccube->numpanelindices;
         }
			
         curStartRow += destcube->numpointindices;
      }
      
      curStartCol = 0;
      for (localSrc = 0; localSrc < numalllocal; localSrc++) {
         srccube = alllocalcubes[localSrc];
         SMatrix_copypiece(Psub, curStartCol, curStartCol,
                          srccube->D_double, 0, 0,
                          srccube->numpointindices, srccube->numpanelindices);
         curStartCol += srccube->numpanelindices;
      }

      // subtract diagonal
      for (i = 0; i < numlocalpoints; i++) {
         diag[i] = Psub[i][i];
         Psub[i][i] = 0.0;
      }
      
      // scale by areas (here we do row scaling because we haven't
      // transposed the double layer operator yet)
      for (i = 0; i < numlocalpoints; i++)
         SVector_scale(Psub[i], areas[i] / (4.0 * M_PI * idiel), numlocalpoints);

      // transpose
      SMatrix_transpose(&Psub, numlocalpoints, numlocalpanels);

      // put correct stuff on diagonal
      // this will break for cavities, need to handle with care like in solv_ecf...
      for (i = 0; i < numlocalpoints; i++)
         Psub[i][i] = (-odiel / ((odiel - idiel) * idiel) + (diag[i]) / (4.0*M_PI*idiel)) * areas[i];

		//		SMatrix_writefile("Psub.m", Psub, numlocalpoints, numlocalpanels);

      // invert Psub
      SMatrix_pseudoinverse(Psubinv, Psub, numlocalpoints, numlocalpanels);


		
      // get rows of Psub corresponding to panels in me
      for (i = 0; i < cube->numpointindices; i++) {
         curpanel = 0;
         for (j = 0; j < numalllocal; j++) {
            for (curcol = 0; curcol < alllocalcubes[j]->numpanelindices; curcol++) {
               Preconditioner_set(preconditioner, cube->pointindices[i],
                                  alllocalcubes[j]->panelindices[curcol],
                                  (real)Psubinv[i][curpanel++]);
            }
         }
      }
/* 		printf("Indices: "); */
/* 		for (i = 0; i < cube->numpointindices; i++) { */
/* 		  curpanel = 0; */
/* 		  printf("%d     ", cube->pointindices[i]); */
/* 		  for (j = 0; j < numalllocal; j++) { */
/* 			 for (curcol = 0; curcol < alllocalcubes[j]->numpanelindices; curcol++) { */
/* 				printf("%d ", alllocalcubes[j]->panelindices[curcol]); */
/* 			 } */
/* 		  } */
/* 		  printf("\n"); */
/* 		} */
/* 		exit(-1); */
		
      // clean up, go home
      SMatrix_free(Psub);
      SMatrix_free(Psubinv);
      free(alllocalcubes);
      free(diag);
      free(areas);
   } else {
      for (cx = 0; cx <= 1; cx++)
         for (cy = 0; cy <= 1; cy++)
            for (cz = 0; cz <= 1; cz++)
               if (cube->children[cx][cy][cz] != NULL) 
                  Preconditioner_fill_overlap_recurse(preconditioner, cube->children[cx][cy][cz],
                                                      panels, idiel, odiel);
   }
}

// this is stolen out of a bunch of functions in Cube.c
void Preconditioner_do_leaf_leaf_translation(SMatrix mat, Cube cube, Cube srccube) {
  unsigned int i, j;
  ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;
  unsigned int halfdimension = 1 + 2 * LOCAL;
  int dx = cube->indices[0] - srccube->indices[0];
  int dy = cube->indices[1] - srccube->indices[1];
  int dz = cube->indices[2] - srccube->indices[2];
  SVector q = SVector_allocate(srccube->numpanelindices);
  SMatrix matT = SMatrix_allocate(srccube->numpanelindices, cube->numpointindices);
  unsigned int* gridpoints = cube->tree->gridpointsperlevel;
  unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
  unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
  SVector VT_q = SVector_allocate(srccube->rowrank);
  SVector PV_VT_q = SVector_allocate(gp3);
  SVector IFFT_sum_T_FFT_PV_VT_q = SVector_allocate(gp3);
  SVector UTI_IFFT_sum_T_FFT_PV_VT_q = SVector_allocate(cube->columnrank);
  ComplexSVector T = Tprecomputed[cube->level][dx+halfdimension][dy+halfdimension][dz+halfdimension];

  for (j = 0; j < srccube->numpanelindices; j++) {
	 SVector_zero(VT_q, srccube->rowrank);
	 SVector_zero(PV_VT_q, gp3);
	 ComplexSVector_zero(srccube->FFT_PV_VT_q, padgridsize);
	 ComplexSVector_zero(cube->sum_T_FFT_PV_VT_q, padgridsize);
	 SVector_zero(IFFT_sum_T_FFT_PV_VT_q, gp3);
	 SVector_zero(UTI_IFFT_sum_T_FFT_PV_VT_q, cube->columnrank);
	 SVector_zero(q, srccube->numpanelindices);
	 q[j] = 1.0;
	 SMatrix_multiplyvector(VT_q, srccube->VTsrc, q, srccube->rowrank, srccube->numpanelindices);
	 SMatrix_multiplyvector(PV_VT_q, srccube->PV_double, VT_q, gp3, srccube->rowrank);
	 FFT_forwardGridTransform(cube->level, gridpoints[cube->level], PV_VT_q, srccube->FFT_PV_VT_q);
 	 ComplexSVector_addelementmultiplyvector(cube->sum_T_FFT_PV_VT_q, T, srccube->FFT_PV_VT_q, padgridsize);
	 FFT_backwardGridTransform(cube->level, gridpoints[cube->level], cube->sum_T_FFT_PV_VT_q, IFFT_sum_T_FFT_PV_VT_q);
	 SMatrix_multiplyvector(UTI_IFFT_sum_T_FFT_PV_VT_q, cube->UTI, IFFT_sum_T_FFT_PV_VT_q, cube->columnrank, gp3);
	 SMatrix_multiplyvector(matT[j], cube->Udest, UTI_IFFT_sum_T_FFT_PV_VT_q, cube->numpointindices, cube->columnrank);
  }
  SVector_free(IFFT_sum_T_FFT_PV_VT_q);
  SVector_free(UTI_IFFT_sum_T_FFT_PV_VT_q);
  SMatrix_transpose(&matT, srccube->numpanelindices, cube->numpointindices);
  SMatrix_copy(mat, matT, cube->numpointindices, srccube->numpanelindices);
  SVector_free(q);
  SMatrix_free(matT);
}

void Preconditioner_fill_overlap_ecf_qual_cav(Preconditioner preconditioner, Tree tree,
                                              Panel* panels, unsigned int numpoints, unsigned int numpanels,
                                              real idiel, real odiel) {
   if (LOCAL != 1) {
      printf("Preconditioner_fill_overlap_ecf_qual_cav:\n");
      printf("\tneeds to be compiled with LOCAL = 1 for now!\n");
      exit(-1);
   }

   Preconditioner_fill_overlap_recurse(preconditioner, tree->root, panels, idiel, odiel);
}

void Preconditioner_multiply(Vector Px, Preconditioner preconditioner, Vector x, unsigned int numpanels) {
   unsigned int i, j;
   Vector_zero(Px, numpanels);
   for (i = 0; i < numpanels; i++) {
/* 	  printf("%d nonzero elements in row\n", preconditioner->numelements[i]); */
	  for (j = 0; j < preconditioner->numelements[i]; j++) {
		 /* 		  printf("%f  ", preconditioner->P[i][j].value); */
		 Px[preconditioner->P[i][j].row] += preconditioner->P[i][j].value *  x[i];
	  }
/* 		printf("\n"); */
   }
}

