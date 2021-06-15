#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"
#include "LinConstrained.h"
#include "BoxConstrained.h"
#include "Overlap.h"

typedef enum { UNCONSTRAINED, LINCONSTRAINED, BOXCONSTRAINED } problemType;
SIZentry* SIZentries; 
unsigned int numSIZentries;
CRGentry* CRGentries;
unsigned int numCRGentries;
char variablechain = 'V';
char fixedligandchain = 'L';
char fixedreceptorchain = 'R';
unsigned int num_GMRES_iter;
int saveGMRES;
real tol;

int main(int argc, char* argv[]) {
  tol = 1e-5;
  //  saveGMRES = 1;
    usequalocation = 1; // for PNAS paper, no salt

  real penaltyScale, tolerance;
  FILE *INPUTFILE;
   unsigned int i,j;
	unsigned int breakout;
	unsigned int numrows, numcols;
	float junk;
	char *s;
	char line[128];
	char *filename, *filename2, *type, *variable, *value, *end, *varname;
	Vector constraintVectorLoaded = NULL;
	Vector upperBounds = NULL;
	Vector lowerBounds = NULL;
	Matrix constraintMatrixLoaded = NULL;
	Matrix penaltyMatrix = NULL;
	problemType curProblemType = UNCONSTRAINED; // default
   PBEproblem bound, unbound;
   UnconstrainedProblem up;
	LinConstrainedProblem lcp;
	BoxConstrainedProblem bcp;
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

	// always use fast M3 products
	accelerateM3 = 1;  // this is a global var from FFTSVDpbeAPI.c

   bound = PBEproblem_allocate(argv[4], argv[5],0);
   unbound   = PBEproblem_allocate(argv[6], argv[7],0);

	Matrix Lbhat = Matrix_allocate(bound->numvariablecharges, bound->numvariablecharges);
	Matrix Luhat = Matrix_allocate(bound->numvariablecharges, bound->numvariablecharges);

	PBEproblem_generateLhatOnTheFly(bound, Lbhat);
	PBEproblem_generateLhatOnTheFly(unbound, Luhat);
	
	Optimizer_computePreconditionerHessian(bound, unbound, Lbhat);
	//	Matrix_writefile("Lhat.m", Lbhat, bound->numvariablecharges, bound->numvariablecharges);

	Matrix_free(Lbhat);
	Matrix_free(Luhat);
	PBEproblem_free(bound);
	PBEproblem_free(unbound);
	return 0;
}
