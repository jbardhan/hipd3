#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"

// true globals--shouldn't be changing between bound and unbound
// states (otherwise we're effectively changing the parameter
// set... and i'm going to assert that my code will not support such
// activity).
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
   PBEproblem bound, unbound;
   unsigned int n_c, columnCount, i;
   unsigned int startcol, endcol;
   Vector chargeVec, unboundReactPot, boundReactPot;
   Matrix Hessian,ligand,complexL;
	usequalocation = 1;
	accelerateM3 = 1;
   if (argc < 13) {
      printf("Usage:\n\t%s <param> <siz> <boundCRG> <boundPDB> <boundSRF> <unboundPDB> <unboundSRF> <startcol> <endcol> <Hessname> <varchain> <fixedligandchain> <receptorchain>\n", argv[0]);
      exit(-1);
   }
	char* foo = argv[11];
	variablechain = foo[0];
	foo = argv[12];
	fixedligandchain = foo[0];
	foo = argv[13];
	fixedreceptorchain = foo[0];

	printf("Ligand variable chain is %c\n", variablechain);
	printf("Receptor chain is %c\n", fixedreceptorchain);
	printf("Ligand fixed chain is %c\n", fixedligandchain);
	
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
   unbound = PBEproblem_allocate(argv[6], argv[7],1);
	//	QualocationOperator_writematlabfile("A2u.m", unbound->qualocationoperator, unbound->numtotalsurfacevariables);
   n_c = unbound->numvariablecharges;
   printf("numvarcharges = %d\n", n_c);
   bound = PBEproblem_allocate(argv[4], argv[5],1);
	//	QualocationOperator_writematlabfile("A2b.m", bound->qualocationoperator, bound->numtotalsurfacevariables);

   chargeVec = Vector_allocate(n_c);
   Hessian = Matrix_allocate(n_c, n_c);
   unboundReactPot = Vector_allocate(n_c);
   boundReactPot = Vector_allocate(n_c);
   startcol = atoi(argv[8]);
   endcol = atoi(argv[9]);
   printf("doing columns %d to %d\n", startcol, (endcol>n_c)?n_c:endcol);
   unsigned int total_GMRES_iter = 0;
   for (columnCount = startcol; columnCount < ((endcol>=n_c)?n_c:endcol+1);
           columnCount++) { 
      for (i = 0; i < n_c; i++)
         chargeVec[i] = 0.0;
      chargeVec[columnCount] = 1.0;

      PBEproblem_setVariableChargeVector(unbound, chargeVec);
		PBEproblem_solve(unbound);
      PBEproblem_getVariableReactionPotentials(unbound, unboundReactPot);
		total_GMRES_iter += num_GMRES_iter; // not counting the much smaller unbound ones.
		num_GMRES_iter = 200;
		
      PBEproblem_setVariableChargeVector(bound, chargeVec);
      PBEproblem_solve(bound);
      PBEproblem_getVariableReactionPotentials(bound, boundReactPot);

		printf("column %d is done: \n", columnCount);
      for (i = 0; i < n_c; i++) {
		  Hessian[i][columnCount] = boundReactPot[i] - unboundReactPot[i];
		  printf("%f  ", Hessian[i][columnCount]);
      }
      printf("\n");
      Vector_zero(boundReactPot, n_c);
      Vector_zero(unboundReactPot, n_c);
  }
   Matrix_writefile(argv[10], Hessian, n_c, n_c);
   Vector_free(chargeVec);
   Matrix_free(Hessian);
   Vector_free(unboundReactPot);
   Vector_free(boundReactPot);
	printf("Success.  Hessian calculation required %d bound-state dense MV products.\n", total_GMRES_iter);
   return 0;
}
