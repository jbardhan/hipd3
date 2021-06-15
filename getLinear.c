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
   unsigned int i;
   PBEproblem bound, unbound;
   UnconstrainedProblem up;
	usequalocation = 1;
	accelerateM3 = 1;
   Vector optimalCharges;
   Vector RHS;
   if (argc < 12) {
      printf("Usage:\n\t%s <param> <siz> <boundCRG> <boundPDB> <boundSRF> <unboundPDB> <unboundSRF> <linearTerm file> <varchain> <fixedvarchain> <receptorchain>\n", argv[0]);
      exit(-1);
   }
   setlinebuf(stdout);

	char* foo = argv[9];
	variablechain = foo[0];
	foo = argv[10];
	fixedligandchain = foo[0];
	foo = argv[11];
	fixedreceptorchain = foo[0];
	printf("chains are %c  %c  %c\n", variablechain, fixedligandchain, fixedreceptorchain);

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

   bound = PBEproblem_allocate(argv[4], argv[5],1);
   unbound   = PBEproblem_allocate(argv[6], argv[7],1);
   UnconstrainedProblem_setPBEproblems(up, bound, unbound);
   RHS = Vector_allocate(2 * (bound->numtotalpanels + unbound->numtotalpanels)+bound->numvariablecharges);

   UnconstrainedProblem_setupRHS(up, RHS);

	FILE *file = NULL;
	file = fopen(argv[8], "w");
	if (file == NULL) {
	  printf("Error opening linear term file! Dying.\n");
	  return -1;
	}
   for (i = 0; i < bound->numvariablecharges; i++)
      fprintf(file, "%f\n", up->linearTerm[i]);
	fclose(file);
	
   UnconstrainedProblem_free(up);
   PBEproblem_free(bound);
   PBEproblem_free(unbound);

   Vector_free(RHS);
   return 0;
}

