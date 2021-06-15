#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"
#include "LinConstrained.h"
#include "BoxConstrained.h"

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
   BoxConstrainedProblem bcp;
   Matrix A_c;
   Vector b;
   Vector lowerBounds, upperBounds;
   Vector optimalCharges, lagrangeMultipliers;

   Vector testq, testslack;
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

   bcp = BoxConstrainedProblem_allocate();

   // load bound, unbound PBE problems
   accelerateM3 = 1; // so A3 products will be done via FFTSVD accel
   bound = PBEproblem_allocate(argv[4], argv[5]);
   unbound = PBEproblem_allocate(argv[6], argv[7]);
   BoxConstrainedProblem_setPBEproblems(bcp, bound, unbound);

   // initialize constraints
   A_c = Matrix_allocate(1, bound->numvariablecharges);
   for (i = 0; i < bound->numvariablecharges; i++)
      A_c[0][i] = 1.0;
   b = Vector_allocate(1);
   b[0] = -2.0;  
   lowerBounds = Vector_allocate(bound->numvariablecharges);
   upperBounds = Vector_allocate(bound->numvariablecharges);
   for (i = 0; i < bound->numvariablecharges; i++) {
      lowerBounds[i] = -0.85;
      upperBounds[i] = +0.85;
   }
   BoxConstrainedProblem_setConstraints(bcp, 1, A_c, b, lowerBounds, upperBounds);

   // MAIN CODE

   optimalCharges = Vector_allocate(bound->numvariablecharges);
   lagrangeMultipliers = Vector_allocate(1); //hard coded numconstraints = 1
   BoxConstrainedProblem_solve(bcp, optimalCharges, lagrangeMultipliers);

   printf("The optimal charge vector is:\n");
   for (i = 0; i < bound->numvariablecharges; i++)
      printf("%f\n", optimalCharges[i]);

   printf("The Lagrange multiplier vector is:\n");
   for (i = 0; i < bcp->numconstraints; i++)
      printf("%f\n", lagrangeMultipliers[i]);


   // END MAIN CODE
   
   // clean up, go home
   BoxConstrainedProblem_free(bcp);
   PBEproblem_free(bound);
   PBEproblem_free(unbound);
   Vector_free(optimalCharges);
   Vector_free(lagrangeMultipliers);
   Vector_free(lowerBounds);
   Vector_free(upperBounds);
   Matrix_free(A_c);
   Vector_free(b);
}
