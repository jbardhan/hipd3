#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"

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
   UnconstrainedProblem up;
   Vector optimalCharges;
   Vector RHS;
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
   up = UnconstrainedProblem_allocate();
   accelerateM3 = 1; //global var from FFTSVDpbeAPI.c : tells
                     //surfaceoperator to use fast A3 product rather
                     //than doing everythign directly
   bound = PBEproblem_allocate(argv[4], argv[5]);
   unbound   = PBEproblem_allocate(argv[6], argv[7]);
   UnconstrainedProblem_setPBEproblems(up, bound, unbound);
   RHS = Vector_allocate(2 * (bound->numtotalpanels + unbound->numtotalpanels)+bound->numvariablecharges);
   optimalCharges = Vector_allocate(bound->numvariablecharges);

   UnconstrainedProblem_setupRHS(up, RHS);
   printf("The linear term is:\n");
   for (i = 0; i < bound->numvariablecharges; i++)
      printf("%f\n", up->linearTerm[i]);

   UnconstrainedProblem_solve(up, optimalCharges);
   
   printf("The optimal charge vector is:\n");
   for (i = 0; i < bound->numvariablecharges; i++)
      printf("%f\n", optimalCharges[i]);


   UnconstrainedProblem_free(up);
   PBEproblem_free(bound);
   PBEproblem_free(unbound);
   Vector_free(optimalCharges);
   Vector_free(RHS);
   return 0;
}

