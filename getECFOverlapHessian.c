#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Overlap.h"

// this file started out as a copy of FFTSVDsolvecfqualcav
#define CONVERSION 332.112
SIZentry* SIZentries; 
unsigned int numSIZentries;
CRGentry* CRGentries;
unsigned int numCRGentries;
char variablechain = 'V';
char fixedligandchain = 'L';
char fixedreceptorchain = 'R';
int saveGMRES;
real tol;

//#define M1M3

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
   Matrix Psub, Psubinv;
   unsigned int numlocalpoints = 0;
   unsigned int numlocalpanels = 0;
   unsigned int numalllocal;
   Cube *alllocalcubes;
   unsigned int i, j, curpanel;
   unsigned int curStartCol, curStartRow, curcol;
   unsigned int localSrc, localDest;
   real *diag, *areas;
   Cube srccube, destcube;
   
   if (cube->leaf) {
      Preconditioner_fill_overlap_get_adjacent_cubes(cube, &numalllocal, &(alllocalcubes));
      for (i = 0; i < numalllocal; i++) {
         numlocalpoints += alllocalcubes[i]->numpointindices;
         numlocalpanels += alllocalcubes[i]->numpanelindices;
      }
      Psub    = Matrix_allocate(numlocalpoints, numlocalpanels);
      Psubinv = Matrix_allocate(numlocalpoints, numlocalpanels);
      diag = (real *)calloc(numlocalpoints, sizeof(real));
      areas = (real *)calloc(numlocalpanels, sizeof(real));
      curpanel = 0;
      for (i = 0; i < numalllocal; i++) {
         for (j = 0; j < alllocalcubes[i]->numpanelindices; j++)
            areas[curpanel++] = panels[alllocalcubes[i]->panelindices[j]]->area;
      }
      // fill submatrix Psub, of P for all localcubes interacting with all local cubes
      curStartRow = 0;
      for (localDest = 0; localDest < numalllocal; localDest++) {
         destcube = alllocalcubes[localDest];
         curStartCol = 0;
         for (localSrc = 0; localSrc < numalllocal; localSrc++) {
            srccube = alllocalcubes[localSrc];
            curcol = destcube->numpanelindices;
            for (i = 0; i < destcube->numlocalcubes; i++) {
               if (srccube == destcube->localcubes[i]) {
                  Matrix_copypiece(Psub, curStartRow, curStartCol,
                                   destcube->D_double, 0, curcol,
                                   destcube->numpointindices, srccube->numpanelindices);
						continue;
               }
					// if we've gotten this far, destcube and srccube are "well separated" via local = 1;
					for (j = 0; j < destcube->numinteractingcubes; j++) {
					  if (srccube == destcube->interactingcubes[j]) {
						 Matrix nearbyCompressedMat = Matrix_allocate(destcube->numpointindices, srccube->numpanelindices);
						 Preconditioner_do_leaf_leaf_translation(nearbyCompressedMat, destcube, srccube);
						 Matrix_copypiece(Psub, curStartRow, curStartCol,
												nearbyCompressedMat, 0, 0,
												destcube->numpointindices, srccube->numpanelindices);
						 Matrix_free(nearbyCompressedMat);
					  }
					}
               curcol += destcube->localcubes[i]->numpanelindices;
            }
            curStartCol += srccube->numpanelindices;
         }
         curStartRow += destcube->numpointindices;
      }
      
      curStartCol = 0;
      for (localSrc = 0; localSrc < numalllocal; localSrc++) {
         srccube = alllocalcubes[localSrc];
         Matrix_copypiece(Psub, curStartCol, curStartCol,
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
         Vector_scale(Psub[i], areas[i] / (4.0 * M_PI * idiel), numlocalpoints);

      // transpose
      Matrix_transpose(&Psub, numlocalpoints, numlocalpanels);

      // put correct stuff on diagonal
      // this will break for cavities, need to handle with care like in solv_ecf...
      for (i = 0; i < numlocalpoints; i++)
         Psub[i][i] = (-odiel / ((odiel - idiel) * idiel) + (diag[i]) / (4.0*M_PI*idiel)) * areas[i];

      // invert Psub
      Matrix_pseudoinverse(Psubinv, Psub, numlocalpoints, numlocalpanels);
      
      // get rows of Psub corresponding to panels in me
      for (i = 0; i < cube->numpointindices; i++) {
         curpanel = 0;
         for (j = 0; j < numalllocal; j++) {
            for (curcol = 0; curcol < alllocalcubes[j]->numpanelindices; curcol++) {
               Preconditioner_set(preconditioner, cube->pointindices[i],
                                  alllocalcubes[j]->panelindices[curcol],
                                  Psubinv[i][curpanel++]);
            }
         }
      }
      // clean up, go home
      Matrix_free(Psub);
      Matrix_free(Psubinv);
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
void Preconditioner_do_leaf_leaf_translation(Matrix mat, Cube cube, Cube srccube) {
  unsigned int i, j;
  ComplexSVector**** Tprecomputed = cube->tree->Tprecomputed;
  unsigned int halfdimension = 1 + 2 * LOCAL;
  int dx = cube->indices[0] - srccube->indices[0];
  int dy = cube->indices[1] - srccube->indices[1];
  int dz = cube->indices[2] - srccube->indices[2];
  Vector q = Vector_allocate(srccube->numpanelindices);
  Matrix matT = Matrix_allocate(srccube->numpanelindices, cube->numpointindices);
  unsigned int* gridpoints = cube->tree->gridpointsperlevel;
  unsigned int padgridsize = (2*gridpoints[cube->level]-1)*(2*gridpoints[cube->level]-1)*((2*gridpoints[cube->level]-1)/2+1);
  unsigned int gp3 = gridpoints[cube->level]*gridpoints[cube->level]*gridpoints[cube->level];
  Vector VT_q = Vector_allocate(srccube->rowrank);
  Vector PV_VT_q = Vector_allocate(gp3);
  Vector IFFT_sum_T_FFT_PV_VT_q = Vector_allocate(gp3);
  Vector UTI_IFFT_sum_T_FFT_PV_VT_q = Vector_allocate(cube->columnrank);
  ComplexSVector T = Tprecomputed[cube->level][dx+halfdimension][dy+halfdimension][dz+halfdimension];

  for (j = 0; j < srccube->numpanelindices; j++) {
	 Vector_zero(q, srccube->numpanelindices);
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
  Matrix_transpose(&matT, srccube->numpanelindices, cube->numpointindices);
  Matrix_copy(mat, matT, cube->numpointindices, srccube->numpanelindices);
  Vector_free(q);
  Matrix_free(matT);
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
      for (j = 0; j < preconditioner->numelements[i]; j++) {
         Px[preconditioner->P[i][j].row] += preconditioner->P[i][j].value *  x[i];
      }
   }
}

int main(int argc, char* argv[]) {
   FILE* vertfile = NULL;
   FILE* facefile = NULL;
   FILE* chargefile = NULL;
   VertFace vertface;
   VertFace* cavvertface;
   Charge charge;
   real svdtol, gmrestol;
   unsigned int gridpoints;
   unsigned int maxpanelsperfinestcube;
   Panel* panels;
   unsigned int numdielpanels;
   unsigned int numcavities = 0;
   unsigned int numpanels;
   Vector3D* centroids;
   unsigned int i, c, currentColumn;
   Vector rhs, q;
   Vector sol;
   Tree tree;
   Matrix Hessian;
   Preconditioner preconditioner;
   Vector phi;
   real energy = 0.0, idiel, odiel;
#ifdef M1M3
   Tree m1m3;
#endif

   if (argc < 11) {
      printf("Usage: %s [surface.vert] [surface.face] [cavitybase] [chargefile] [svd tol] [gmres tol] [gridpoints] [maxpanels] [idiel] [odiel]\n", argv[0]);
      return -1;
   }

   vertface = VertFace_allocate();

   vertfile = fopen(argv[1], "r");

   if (vertfile == NULL) {
      perror("Error opening vertices file");
      return -2;
   }

   VertFace_readvert(vertface, vertfile);

   fclose(vertfile);

   facefile = fopen(argv[2], "r");
    
   if (facefile == NULL) {
      perror("Error opening face file");
      return -2;
   }
    
   VertFace_readface_flip(vertface, facefile);

   for (i = 1; ; i++) {
      char filename[256];
      sprintf(filename, "%s_%u.face", argv[3], i);
      if (!access(filename, R_OK))
         numcavities++;
      else
         break;
   }

   printf("NUMCAVITIES: %u\n", numcavities);

   cavvertface = (VertFace*)calloc(numcavities, sizeof(VertFace));

   for (i = 0; i < numcavities; i++) {
      char filename[256];
      sprintf(filename, "%s_%u.vert", argv[3], i+1);
      cavvertface[i] = VertFace_allocate();
      FILE* cavvertfile = fopen(filename, "r");
      VertFace_readvert(cavvertface[i], cavvertfile);
      fclose(cavvertfile);
      sprintf(filename, "%s_%u.face", argv[3], i+1);
      FILE* cavfacefile = fopen(filename, "r");
      VertFace_readface_flip(cavvertface[i], cavfacefile);
      fclose(cavfacefile);
   }

   charge = Charge_allocate();
   
   chargefile = fopen(argv[4], "r");
   
   if (chargefile == NULL) {
      perror("Error opening charge file");
      return -2;
   }

   Charge_read(charge, chargefile);
   
   fclose(chargefile);

   svdtol = atof(argv[5]);
   gmrestol = atof(argv[6]);
   gridpoints = atoi(argv[7]);
   maxpanelsperfinestcube = atoi(argv[8]);
   idiel = atoi(argv[9]);
   odiel = atoi(argv[10]);

   if (odiel < 1.1) {
      VertFace_fix(vertface, 1);
      for (i = 0; i < numcavities; i++)
         VertFace_fix(cavvertface[i], 1);
   }
   else {
      VertFace_fix(vertface, 0);
      for (i = 0; i < numcavities; i++)
         VertFace_fix(cavvertface[i], 0);
   }

   numdielpanels = vertface->numfaces;
   numpanels = numdielpanels;
   for (i = 0; i < numcavities; i++)
      numpanels += cavvertface[i]->numfaces;

   printf("NUMPANELS: %u\n", numpanels);

   panels = (Panel*)calloc(numpanels, sizeof(Panel));

   printf("Constructing and allocating panels... ");
   fflush(stdout);

   VertFace_getpanels(vertface, panels);
   c = 0;
   for (i = 0; i < numcavities; i++) {
      VertFace_getpanels(cavvertface[i], panels + numdielpanels + c);
      c += cavvertface[i]->numfaces;
   }

   printf("done.\n");

   real area = 0.0;
   for (i = 0; i < numpanels; i++)
      area += panels[i]->area;

   printf("Surface Area: %f\n", area);

   centroids = (Vector3D*)calloc(numpanels, sizeof(Vector3D));
   
   for (i = 0; i < numpanels; i++) {
      centroids[i] = Vector3D_allocate();
      Vector3D_copy(centroids[i], panels[i]->centroid);
   }

   rhs = Vector_allocate(numpanels);

   printf("doing %d charges\n", charge->numcharges);
#ifdef M1M3
   m1m3 = Tree_allocate(panels, numpanels, charge->points, charge->numcharges, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, SINGLE_AND_DOUBLE_LAYER_INT, 0.0);
   Tree_lists(m1m3);
   Tree_fill(m1m3);
   Tree_memory(m1m3);
   q = Vector_allocate(charge->numcharges);
#endif

   sol = Vector_allocate(numpanels);

   printf("Constructing and allocating tree... ");
   fflush(stdout);

   tree = Tree_allocate(panels, numpanels, centroids, numpanels, maxpanelsperfinestcube, POISSON_KERNEL, NULL, gridpoints, svdtol, DOUBLE_LAYER_INT, 0.0);

   printf("done.\n");

   printf("Determining local and interacting lists... ");
   fflush(stdout);

   Tree_lists(tree);

   printf("done.\n");

   printf("Filling tree structure... ");
   fflush(stdout);

   Tree_fill(tree);

   printf("done.\n");

   Tree_memory(tree);

   printf("Constructing preconditioner... ");
   fflush(stdout);
   
   preconditioner = Preconditioner_allocate(numpanels, numpanels);

   Preconditioner_fill_overlap_ecf_qual_cav(preconditioner, tree, panels, numpanels,
                                            numpanels, idiel, odiel);
   // do not factor the overlap preconditioner

/*    Preconditioner_fill_diagonal_solv_ecf_qual_cav(preconditioner, panels, numpanels, numpanels, idiel, odiel); */
/*    Preconditioner_factor(preconditioner); */

   printf("done.\n");

   printf("Memory use for P: %u\n", Preconditioner_memory(preconditioner));

   phi = Vector_allocate(charge->numcharges);
   
   Hessian = Matrix_allocate(charge->numcharges, charge->numcharges);
   for (currentColumn = 0; currentColumn < charge->numcharges; currentColumn++) { 
      for (i = 0; i < charge->numcharges; i++)
         charge->charges[i] = 0.0;
      charge->charges[currentColumn] = 1.0;

#ifdef M1M3
      for (i = 0; i < charge->numcharges; i++)
         q[i] = -charge->charges[i] / (4.0 * M_PI * idiel);
      Tree_multiply_transpose(rhs, m1m3, q, DOUBLE_LAYER_INT);
#else
      Charge_makerhs_ecf_qual(rhs, charge, panels, numpanels, idiel, odiel);
#endif

      
      // instead of GMRES here, just apply precond
      //      Preconditioner_solve(sol, preconditioner, rhs);
      Preconditioner_multiply(sol, preconditioner, rhs, numpanels);

#ifdef M1M3
   Tree_multiply(phi, m1m3, sol, SINGLE_LAYER_INT);
   Vector_scale(phi, CONVERSION / idiel, charge->numcharges);
#else
   for (i = 0; i < numpanels; i++) {
      for (c = 0; c < charge->numcharges; c++) {
         real slp, dlp;
         
         slp = Integration(charge->points[c], panels[i], POISSON_KERNEL, NULL, SINGLE_LAYER_INT);
         
         phi[c] += CONVERSION * (sol[i] * slp) / idiel;
      }
   }
#endif

   for (i = 0; i < charge->numcharges; i++)
      Hessian[i][currentColumn] = phi[i];
   }

   Matrix_writefile(argv[11], Hessian, charge->numcharges, charge->numcharges);
#ifdef M1M3
   Tree_free(m1m3);
   Vector_free(q);
#endif
   Vector_free(sol);
   Vector_free(rhs);
   Vector_free(phi);
   for (i = 0; i < numpanels; i++)
      Vector3D_free(centroids[i]);
   free(centroids);
   for (i = 0; i < numpanels; i++)
      Panel_free(panels[i]);
   free(panels);
   Charge_free(charge);
   VertFace_free(vertface);
   for (i = 0; i < numcavities; i++)
      VertFace_free(cavvertface[i]);
   free(cavvertface);

   return 0;

}
