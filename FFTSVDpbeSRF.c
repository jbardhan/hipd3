#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"

// this file is a ripoff of getHessian, only all it does is compute a
// single solvation energy.

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
unsigned int num_GMRES_iter;
int saveGMRES = 0;
real tol;

void doAllCharges(PBEproblem problem) {
  problem->numvariablecharges = problem->numpdbentries;
  problem->numfixedligandcharges = 0;
  problem->numfixedreceptorcharges = 0;
  free(problem->variablechargeindextoglobalindex);
  free(problem->fixedligandchargeindextoglobalindex);
  free(problem->fixedreceptorchargeindextoglobalindex);
  problem->variablechargeindextoglobalindex = (unsigned int *)calloc(problem->numvariablecharges, sizeof(unsigned int));
  unsigned int i;
  for (i =0; i < problem->numvariablecharges; i++) {
	 problem->variablechargeindextoglobalindex[i] = i;
  }
}

int main(int argc, char* argv[]) {
  PBEproblem problem;
  unsigned int i;
  real energy = 0.0;
   Matrix Hessian;
   if (argc < 7) {
      printf("Usage:\n\t%s <param> <siz> <boundCRG> <boundPDB> <boundSRF> doCoulomb? chainsToCharge\n", argv[0]);
      exit(-1);
   }
   setlinebuf(stdout);
   usequalocation = 0;
   // initialization of Tidor Lab data structures
   printf("Reading parameters from %s\n", argv[1]);
   readParams(argv[1]);
   printf("Checking for valid parameters\n");
   checkParams();
   printf("Reading radii from %s\n", argv[2]);
   readSIZ(argv[2], &numSIZentries, &SIZentries);
   printf("Reading charges from %s\n", argv[3]);
   readCRG(argv[3], &numCRGentries, &CRGentries);
	
   num_GMRES_iter = 1000;
   // load bound, unbound PBE problems
	accelerateM3 = 1;
   problem = PBEproblem_allocate(argv[4], argv[5], 1);
	doAllCharges(problem);

	//	QualocationOperator_writematlabfile("A_qual.m", problem->qualocationoperator, problem->numtotalpanels);
	// 	SurfaceOperator_writematlabfile("A_green.m", problem->pbesurfaceoperator, problem->numtotalpanels);

	//	exit(-1);
	
   for (i = 0; i < problem->numpdbentries; i++) {
      problem->globalCharges[i] = problem->pdbentries[i].charge;
	}
	
	PBEproblem_solve(problem);

/*    printf("reaction potentials = \n"); */
	unsigned int j;
	for (i = 0; i < problem->numpdbentries; i++)
	  energy += problem->globalPhiReact[i] * problem->globalCharges[i] * .5; // .592 already accounted for
	printf("total solvation energy = %f\n", energy);

	real CoulombEnergy = 0.0;
	Vector coulombPotential = Vector_allocate(problem->numpdbentries);
	Vector3D* pointList = problem->qualocationoperator->charges->points;
	Vector charges = problem->qualocationoperator->charges->charges;
	if (argc > 5) {
	  for (i = 0; i < problem->numpdbentries; i++) {
		 //		 printf("i = %d, ce = %f\n", i, CoulombEnergy);
		 coulombPotential[i] = 0.0;
		 for (j = 0; j < problem->numpdbentries; j++) {
			if (i != j) {
			  coulombPotential[i] += problem->globalCharges[j] / (Vector3D_distance(pointList[i], pointList[j]));
			}
		 }
		 coulombPotential[i] *= 0.592 * KT_CONVERSION / innerdielectric;
/* 		 printf("coulombPotential[i] = %f  q[i] = %f\n", coulombPotential[i], */
/* 				  problem->globalCharges[i]); */
		 CoulombEnergy += 0.5 * coulombPotential[i] *problem->globalCharges[i];
	  }
	  //	  printf("i = %d, ce = %f\n", i, CoulombEnergy);
	}
	Vector_free(coulombPotential);
	printf("total Coulomb energy = %f\n", CoulombEnergy);

   PBEproblem_free(problem);

	return 0;
}
