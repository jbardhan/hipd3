#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"
#include "LinConstrained.h"

SIZentry* SIZentries; 
unsigned int numSIZentries;
CRGentry* CRGentries;
unsigned int numCRGentries;
char variablechain = 'V';
char fixedligandchain = 'L';
char fixedreceptorchain = 'R';
unsigned int num_GMRES_iter = 1000;
int saveGMRES;
real tol;

int main(int argc, char* argv[]) {
  tol = 1e-5;
  saveGMRES = 0;
  usequalocation = 1;
   unsigned int i;
   PBEproblem bound, unbound;
   LinConstrainedProblem lcp;
   Matrix A_c;
   Vector b;
   Vector optimalCharges, lagrangeMultipliers;
   
   if (argc < 9) {
      printf("Usage:\n\t%s <param> <siz> <boundCRG> <boundPDB> <boundSRF> <unboundPDB> <unboundSRF> <optionsFile>\n", argv[0]);
      exit(-1);
   }
   setlinebuf(stdout);
   
   // initialization of Tidor Lab data structures
   printf("Reading parameters from %s\n", argv[1]);
   readParams(argv[1]);
   printf("Checking for valid parameters\n");
   checkParams();
   printf("Reading radii from %s\n", argv[2]);
   readSIZ(argv[2], &numSIZentries, &SIZentries);
   printf("Reading charges from %s\n", argv[3]);
   readCRG(argv[3], &numCRGentries, &CRGentries);
   
   
   // load bound, unbound PBE problems
   accelerateM3 = 1; // so A3 product will be done via FFTSVD accel
   lcp = LinConstrainedProblem_allocate();
   bound = PBEproblem_allocate(argv[4], argv[5]);
   printf("done allocating bound problem\n");
   unbound   = PBEproblem_allocate(argv[6], argv[7]);
   printf("done allocating unbound problem\n");
   LinConstrainedProblem_setPBEproblems(lcp, bound, unbound);
   
   A_c = Matrix_allocate(1, lcp->bound->numvariablecharges);
   for (i =0; i < bound->numvariablecharges; i++)
      A_c[0][i] = 1.0;
   b = Vector_allocate(1);
   b[0] = -2.0;  
   lagrangeMultipliers = Vector_allocate(1);
   LinConstrainedProblem_setConstraints(lcp, 1, A_c, b);
/*    LinConstrainedProblem_loadConstraints(argv[8]); */

   optimalCharges = Vector_allocate(bound->numvariablecharges);
   LinConstrainedProblem_solve(lcp, optimalCharges, lagrangeMultipliers);

   printf("The optimal charge vector is:\n");
   for (i = 0; i < bound->numvariablecharges; i++)
      printf("%f\n", optimalCharges[i]);

   printf("The Lagrange multiplier vector is:\n");
   for (i = 0; i < lcp->numconstraints; i++)
      printf("%f\n", lagrangeMultipliers[i]);
   
   Matrix_free(A_c);
   Vector_free(b);
   Vector_free(lagrangeMultipliers);   

   LinConstrainedProblem_free(lcp);
   PBEproblem_free(bound);
   PBEproblem_free(unbound);
   Vector_free(optimalCharges);
   return 0;
}

