#ifndef __UNCONSTRAINED_H__
#define __UNCONSTRAINED_H__

#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"

typedef struct _UnconstrainedProblem {
   PBEproblem bound;
   PBEproblem unbound;

   Vector linearTerm;
   Matrix penaltyMatrix;
   real penalty;
   Matrix LhatInv; // for preconditioner
   Matrix Lhat;

   Matrix leftSingularVectors;
   Matrix rightSingularVectors;
   Vector singularValues;
  
  penaltyToleranceType penaltyType;
  real maxEigenvalue; // this is approximate!
  real tolerance; 
  real looserTol;
  real tighterTol;
} _UnconstrainedProblem;

typedef struct _UnconstrainedProblem* UnconstrainedProblem;

UnconstrainedProblem UnconstrainedProblem_allocate();
void UnconstrainedProblem_free(UnconstrainedProblem up);
void UnconstrainedProblem_setPBEproblems(UnconstrainedProblem up, PBEproblem bound, PBEproblem unbound);
void UnconstrainedProblem_setPenalty(UnconstrainedProblem up, penaltyToleranceType penaltyType,
												 real maxEig, real penalty, Matrix penaltyMatrix,
												 Matrix leftSingularVectors, Matrix rightSingularVectors,
												 Vector singularValues,
												 real tolerance, real looserTol, real tighterTol);
void UnconstrainedProblem_solve(UnconstrainedProblem up, Vector optimalCharges,
										  Vector** approxCharges, int *numAdjustments, int* lowerAdjustment,
										  int* higherAdjustment);
void UnconstrainedProblem_getApproxSWMSoln(UnconstrainedProblem up, Vector newsoln, Vector oldsoln,
														 unsigned int orig, unsigned int newInd);
void UnconstrainedProblem_GMRES(UnconstrainedProblem up, Vector soln, Vector RHS);
void UnconstrainedProblem_setupPreconditioner(UnconstrainedProblem up);
void UnconstrainedProblem_setupRHS(UnconstrainedProblem up, Vector RHS);
void UnconstrainedProblem_operatorMultiply(UnconstrainedProblem up, Vector Ax, Vector x);
void UnconstrainedProblem_preconditionerMultiply(UnconstrainedProblem up, Vector Px, Vector x);
void UnconstrainedProblem_preconditionerMultiplyOrig(UnconstrainedProblem up, Vector Px, Vector x);
void UnconstrainedProblem_preconditionerMultiplyLowerTriang(UnconstrainedProblem up, Vector Px, Vector x);
void UnconstrainedProblem_loadNewChargeDistribution(UnconstrainedProblem up, char *PDBfilename, char *CRGfilename);

void UnconstrainedProblem_operatorSave(UnconstrainedProblem up, char *filename);
void UnconstrainedProblem_preconditionerSave(UnconstrainedProblem up, char *filename);

void UnconstrainedProblem_join(UnconstrainedProblem up,
                        Vector dest, Vector q_temp, Vector phi_b_temp, Vector phiu_temp);
void UnconstrainedProblem_split(UnconstrainedProblem up,
                         Vector q_temp, Vector phi_b_temp, Vector phiu_temp, Vector src);


#endif 
