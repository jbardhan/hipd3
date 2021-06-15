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
unsigned int num_GMRES_iter = 1000;
int saveGMRES;
real tol;

int main(int argc, char* argv[]) {
  penaltyToleranceType penaltyType = RELATIVE_RAYLEIGH_TOL;
  real maxEig;
  tol = 1e-5;
  //  saveGMRES = 1;
  usequalocation = 1; // for PNAS paper, no salt

  real penaltyScale, tolerance;
  real looserTol, tighterTol;
  Vector *approxCharges;
  int numCharges;
  int *numAddedLooser;
  int *numAddedTighter;
  CRGentries = NULL;
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
	Matrix rightSingularVectors = NULL;
	Matrix leftSingularVectors = NULL;
	Vector singularValues = NULL;
	problemType curProblemType = UNCONSTRAINED; // default
   PBEproblem bound, unbound;
   UnconstrainedProblem up;
	LinConstrainedProblem lcp;
	BoxConstrainedProblem bcp;
	
   Vector optimalCharges;
	Vector lagrangeMultipliers;
   Vector RHS_up, RHS_lp, RHS_bp;
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
   unbound   = PBEproblem_allocate(argv[6], argv[7],1);
   bound = PBEproblem_allocate(argv[4], argv[5],1);

   optimalCharges = Vector_allocate(bound->numvariablecharges);

	// set up unconstrained problem
   up = UnconstrainedProblem_allocate();
   UnconstrainedProblem_setPBEproblems(up, bound, unbound);

	// set up linconstrained problem
	lcp = LinConstrainedProblem_allocate();
	LinConstrainedProblem_setPBEproblems(lcp, bound, unbound);

	// set up boxconstrained problem
	bcp = BoxConstrainedProblem_allocate();
	BoxConstrainedProblem_setPBEproblems(bcp, bound, unbound);

	breakout = 0;
	do {
	  printf("HIPD> ");
	  // get line of input
	  fgets(line, sizeof(line), stdin);

	  // parse input
	  for (s = line; (s = strtok(s, " \t\n")) != NULL; s = NULL) {
		 if (s[0] == '#') {
			printf("Info: Skipping comment character and the remainder of the line\n");
			break;
		 } else if (strcmp(s, "load")==0) {/** LOAD **********************************************/
			// options: eq constraint, ineq, distribution
			type = strtok(NULL, " \t\n");
			if (type == NULL) {
			  printf("ERROR: load command needs a type.\n");
			}
			filename = strtok(NULL, " \t\n");
			if (filename == NULL) {
			  printf("ERROR: load command needs a filename.\n");
			  break;
			}
			INPUTFILE = fopen(filename, "r");
			if (INPUTFILE == NULL) {
			  printf("ERROR: Could not open input file!\n");
			  break;
			}
			printf("loading %s...\n", filename);
			if (strcmp(type, "eq") == 0) {  // WE ARE NOT DOING ANY ERROR OR EOF CHECKING!
			  fscanf(INPUTFILE, "%d %d", &numrows, &numcols);
			  printf("%d by %d constraint matrix\n", numrows, numcols);
			  if (constraintMatrixLoaded != NULL) {
				 Matrix_free(constraintMatrixLoaded);
				 Vector_free(lagrangeMultipliers);
			  }
			  lagrangeMultipliers = Vector_allocate(numrows);
			  constraintMatrixLoaded = Matrix_allocate(numrows, numcols);
			  for (i = 0; i < numrows; i++) {
				 for (j = 0; j < numcols; j++) {
					fscanf(INPUTFILE, "%f ", &junk);
					constraintMatrixLoaded[i][j] = (real) junk;
					printf("%f  ", constraintMatrixLoaded[i][j]);
				 }
			  }
			  printf("Info: Constraint matrix successfully loaded.\n");
			  if (constraintVectorLoaded != NULL){
				 Vector_free(constraintVectorLoaded);
			  }
			  constraintVectorLoaded = Vector_allocate(numrows);
			  for (i = 0; i < numrows; i++) {
				 fscanf(INPUTFILE, "%f ", &junk);
				 constraintVectorLoaded[i] = (real) junk;				 
				 printf("%f  ", constraintVectorLoaded[i]);
			  }
			  LinConstrainedProblem_setConstraints(lcp, numrows, constraintMatrixLoaded, constraintVectorLoaded);
			  if (lowerBounds != NULL) {
				 BoxConstrainedProblem_setConstraints(bcp, numrows, constraintMatrixLoaded,
																  constraintVectorLoaded, lowerBounds, upperBounds);
			  }
			  printf("Info: Constraint vector successfully loaded.\n");
			}
			if (strcmp(type, "ineq") ==0) {
			  lowerBounds = Vector_allocate(bound->numvariablecharges);
			  upperBounds = Vector_allocate(bound->numvariablecharges);
			  for (i = 0; i < bound->numvariablecharges; i++) {
				 fscanf(INPUTFILE, "%f ", &junk);
				 lowerBounds[i] = (real) junk;
			  }
			  for (i = 0; i < bound->numvariablecharges; i++) {
				 fscanf(INPUTFILE, "%f ", &junk);
				 upperBounds[i] = (real) junk;
			  }
			  BoxConstrainedProblem_setConstraints(bcp, numrows, constraintMatrixLoaded,
																constraintVectorLoaded,
																lowerBounds, upperBounds);
			}

			// new distribution -- this doesn't work as of 2/26/07
			if (strcmp(type, "dist") == 0) {
			  filename2 = strtok(NULL, " \t\n");
			  if (filename2 == NULL) {
				 printf("ERROR: load dist command needs TWO filenames\n");
			  }
			  else {
				 PBEproblem_loadNewChargeDistribution(bound, filename, filename2);
				 PBEproblem_loadNewChargeDistribution(unbound, filename, filename2);
				 Optimizer_updateChargeDistributions(bound, unbound);
			  }
			}

			// some quadratic penalty matrix
			if (strcmp(type,"penaltyMatrix") == 0) {
			  if (penaltyMatrix != NULL) {
				 Matrix_free(penaltyMatrix);
			  }
			  
			  penaltyMatrix = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  leftSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  rightSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  singularValues = Vector_allocate(up->bound->numvariablecharges);
			  for (i = 0; i < up->bound->numvariablecharges; i++) {
				 for (j = 0; j < up->bound->numvariablecharges; j++) {
					fscanf(INPUTFILE, "%f ", &junk);
					penaltyMatrix[i][j] = (real) junk;
				 }
			  }
			  printf("this feature is currently de-implemented!\n");
/* 			  Unconstrainedproblem_setPenalty(up, penaltyMatrix); */
/* 			  LinConstrainedProblem_setPenalty(lcp, penaltyMatrix); */
/* 			  BoxConstrainedProblem_setPenalty(bcp, penaltyMatrix); */
			} 

			fclose(INPUTFILE);
			/* END LOAD */

		 } else if (strcmp(s, "save")==0) {/* SAVE ***********************************************/
			// options: charges, multipliers, GMRES basis + Hessenberg
			type = strtok(NULL, " \t\n");
			if (type == NULL) {
			  printf("ERROR: save command needs a type.\n");
			}

			filename = strtok(NULL, " \t\n");
			if (filename == NULL) {
			  printf("ERROR: save command needs a filename\n");
			  break;
			}
			if (strcmp(type, "operator") == 0) {
			  UnconstrainedProblem_operatorSave(up, filename);
			  break;
			} else if (strcmp(type, "preconditioner") == 0 ) {
			  UnconstrainedProblem_preconditionerSave(up, filename);
			  break;
			}
			FILE *OUTPUTFILE = fopen(filename, "w");
			varname = strtok(NULL, " \t\n");
			printf("saving %s...\n", filename);
			if (strcmp(type,"optimalq") ==0) {
			  if (varname == NULL) {
				 fprintf(OUTPUTFILE, "optimalq = [");
			  } else {
				 fprintf(OUTPUTFILE, "%s = [", varname);
			  }
			  for (i = 0; i < bound->numvariablecharges; i++) {
				 fprintf(OUTPUTFILE, "%f\n", optimalCharges[i]);
			  }
			} else if (strcmp(type, "eq") == 0 ) {
			  fprintf(OUTPUTFILE, "lagrangeMult = [");
			  for (i = 0; i < numrows ; i++) {
				 fprintf(OUTPUTFILE, "%f\n", lagrangeMultipliers[i]);
			  }
			} else if (strcmp(type, "L") == 0) {
			  Matrix L = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  Optimizer_computeHessian(up->bound, up->unbound, L, 0, up->bound->numvariablecharges + 1);
			  fprintf(OUTPUTFILE, "L = [");
			  for (i = 0; i < up->bound->numvariablecharges; i++) {
				 for (j = 0; j < up->bound->numvariablecharges; j++) {
					fprintf(OUTPUTFILE, "%f ", L[i][j]);
				 }
				 fprintf(OUTPUTFILE, "\n");
			  }
			} else if (strcmp(type, "c") == 0) {
			  if (up->linearTerm == NULL) {
				 up->linearTerm = Vector_allocate(up->bound->numvariablecharges);
				 Optimizer_computeLinearTerm(up->bound, up->unbound, up->linearTerm);
			  }
			  fprintf(OUTPUTFILE, "c = [");
			  for (i = 0; i < up->bound->numvariablecharges; i++) {
				 fprintf(OUTPUTFILE, "%f\n", up->linearTerm[i]);
			  }
			} else if ((strcmp(type, "Ldhat") == 0)||(strcmp(type, "Lohat")==0)) {
			  fprintf(OUTPUTFILE, "%s = [", type);
			  if (strcmp(type,"Ldhat") ==0) {
				 up->bound->useOverlap = 0;
				 up->unbound->useOverlap = 0;
			  }
			  if (strcmp(type,"Lohat") ==0) {
				 up->bound->useOverlap = 1;
				 up->unbound->useOverlap = 1;
			  }
			  Matrix Lhat = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  Optimizer_computePreconditionerHessian(up->bound, up->unbound, Lhat);
			  
			  for (i = 0; i < up->bound->numvariablecharges; i++) {
				 for (j = 0; j < up->bound->numvariablecharges; j++) {
					fprintf(OUTPUTFILE,"%f ", Lhat[i][j]);
				 }
				 fprintf(OUTPUTFILE, "\n");
			  }
			  Matrix_free(Lhat);
			} 
			fprintf(OUTPUTFILE, "];\n");
			fclose(OUTPUTFILE);
			break;
			/* END SAVE */
		 } else if (strcmp(s, "solve")==0) {/* SOLVE ***********************************************/
			// options: none
			num_GMRES_iter = 0;
			printf("solving...\n");
			int numApproxSolns = 0;
			int indexLooser;
			int indexTighter;
			
			if (curProblemType == UNCONSTRAINED) { // really want there to be one solve for linear term only.
			  printf("solving unconstrained problem\n");
			  UnconstrainedProblem_solve(up, optimalCharges,
												  &approxCharges, &numApproxSolns, &indexLooser, &indexTighter);
			  char filenameApprox[100];
			  printf("from %d to %d vectors penalized\n", indexLooser, indexTighter);
			  for (i = 0; i < numApproxSolns; i++) {
				 sprintf(filenameApprox, "unc_%d.m", bound->numvariablecharges - (indexLooser+i));
				 FILE *APPROX = fopen(filenameApprox, "w");
				 fprintf(APPROX, "unc_approx_%d = [", bound->numvariablecharges - (indexLooser+i));
				 for (j = 0; j < bound->numvariablecharges; j++) {
					fprintf(APPROX, "%f\n", (approxCharges[i])[j]);
				 }
				 fprintf(APPROX, "];\n");
				 fclose(APPROX);
			  }
			} else if (curProblemType == LINCONSTRAINED) {
			  // should catch no constraints set
			  printf("solving linconstrained problem\n");
			  LinConstrainedProblem_solve(lcp, optimalCharges, lagrangeMultipliers,
													&approxCharges, &numApproxSolns, &indexLooser, &indexTighter);
			} else if (curProblemType == BOXCONSTRAINED) {
			  // should catch no constraints set
			  printf("solving boxconstrained problem\n");
			  BoxConstrainedProblem_solve(bcp, optimalCharges, lagrangeMultipliers);
			} else {
			  printf("Illegal problem type!\n");
			}
			printf("Problem took %d GMRES iterations to solve.\n", num_GMRES_iter);
			break;
			/* END SOLVE */
		 } else if (strcmp(s, "set")==0) {/* SET ***********************************************/
			// options: GMREStol, saveGMRESbasis, diagonalPenalty, problemType
			variable = strtok(NULL, " \t\n");
			value = strtok(NULL, " \t\n");
			if ((variable == NULL) || (value == NULL)) {
			  printf("ERROR: set command needs a variable.\n");
			}
			if (strcmp(variable,"GMREStol") == 0) {
			  tol = (real)strtof(value, &end);
			} else if (strcmp(variable, "saveGMRESbasis") == 0) {
			  saveGMRES = (int)atoi(value);
			} else if (strcmp(variable, "penalty") == 0) {
			  variable = strtok(NULL, " \t\n"); // this is the tolerance
			  if (value == NULL) {
				 printf("WARNING: penaltyScale taken to be default 10.\n");
				 penaltyScale = 10.;
			  } else {
				 penaltyScale = (real) strtod(value,NULL);
			  }
			  if (variable == NULL) {
				 printf("WARNING: tolerance taken to be default 1e-4.\n");
				 tolerance = 1e-4;
			  } else {
				 tolerance = (real) strtod(variable, NULL);
			  }
			  looserTol = tolerance;
			  tighterTol = tolerance;
/* 			  char filename[100]; */
/* 			  sprintf(filename,"data_1/Lhat_%d.m",up->bound->numvariablecharges); */
			  leftSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  rightSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  singularValues = Vector_allocate(up->bound->numvariablecharges);
			  Matrix Lhat = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  Optimizer_computePreconditionerHessian(up->bound, up->unbound, Lhat);
			  penaltyMatrix = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
			  printf("tolerance = %f\n", tolerance);
			  Optimizer_computePenaltyMatrix(up->bound, up->unbound, up->bound->numvariablecharges, Lhat, penaltyType, tolerance, penaltyScale, penaltyMatrix, leftSingularVectors, rightSingularVectors, singularValues, &maxEig);
			  UnconstrainedProblem_setPenalty(up, penaltyType, penaltyScale, maxEig, penaltyMatrix, leftSingularVectors, rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  LinConstrainedProblem_setPenalty(lcp, penaltyScale, penaltyMatrix, leftSingularVectors, rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  BoxConstrainedProblem_setPenalty(bcp, penaltyScale, penaltyMatrix, leftSingularVectors, rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			} else if (strcmp(variable, "problemType") == 0) {
			  if (strcmp(value, "unc") == 0) {
				 curProblemType = UNCONSTRAINED;
			  } else if (strcmp(value, "lin") == 0) {
				 curProblemType = LINCONSTRAINED;
			  } else if (strcmp(value, "box") == 0) {
				 curProblemType = BOXCONSTRAINED;
			  }
			} else if (strcmp(variable, "looserTol") ==  0) {
			  looserTol = (real)strtod(value, NULL);
			  UnconstrainedProblem_setPenalty(up, penaltyType, penaltyScale, maxEig, penaltyMatrix, leftSingularVectors,
														 rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  LinConstrainedProblem_setPenalty(lcp, penaltyScale, penaltyMatrix, leftSingularVectors,
														  rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  BoxConstrainedProblem_setPenalty(bcp, penaltyScale, penaltyMatrix, leftSingularVectors,
														  rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);

			} else if (strcmp(variable, "tighterTol") ==  0) {
			  tighterTol = (real)strtod(value, NULL);
			  UnconstrainedProblem_setPenalty(up, penaltyType, penaltyScale, maxEig, penaltyMatrix, leftSingularVectors,
														 rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  LinConstrainedProblem_setPenalty(lcp, penaltyScale, penaltyMatrix, leftSingularVectors,
														  rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			  BoxConstrainedProblem_setPenalty(bcp, penaltyScale, penaltyMatrix, leftSingularVectors,
														  rightSingularVectors, singularValues,
														 tolerance, looserTol, tighterTol);
			}
			printf("setting...\n");
			break;
			/* END SET */
		 } else if ((strcmp(s, "quit")==0) || (strcmp(s, "q") == 0) || (strcmp(s, "exit") == 0)) {
			printf("Info: Quitting...\n");
			breakout = 1;
		 } else if (strcmp(s, "fix") == 0) {/* FIX ***********************************************/
			printf("ERROR: fix command is not implemented yet!\n");
			break;
		 } else {
			printf("ERROR: Command not understood.\n");
			break;
		 }
	  }
	} while (! breakout);

	return 0;
}
