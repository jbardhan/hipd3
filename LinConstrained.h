#ifndef __LINCONSTRAINED_H__
#define __LINCONSTRAINED_H__

#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"

typedef struct _LinConstrainedProblem {
   PBEproblem bound;
   PBEproblem unbound;

   Vector linearTerm;
   unsigned int numconstraints;
   Matrix A_c;  // these are the equality constraints
   Vector b;  
   
  Matrix penaltyMatrix;
  real penalty;
  Matrix Lhat; 
  Matrix LhatAndAcInv; // for the preconditioner

  Matrix leftSingularVectors;
  Matrix rightSingularVectors;
  Vector singularValues;
  
  real tolerance; 
  real looserTol;
  real tighterTol;

} _LinConstrainedProblem;

typedef struct _LinConstrainedProblem* LinConstrainedProblem;
// public interface
LinConstrainedProblem LinConstrainedProblem_allocate();
void LinConstrainedProblem_free(LinConstrainedProblem lcp);
void LinConstrainedProblem_loadConstraints(LinConstrainedProblem lcp, char* filename);
void LinConstrainedProblem_setConstraints(LinConstrainedProblem lcp, unsigned int numconstraints,
                                          Matrix A_c, Vector b);
void LinConstrainedProblem_setPenalty(LinConstrainedProblem lcp, real penalty, Matrix penaltyMatrix,
												  Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues,
												  real tolerance, real looserTol, real tighterTol);
void LinConstrainedProblem_setPBEproblems(LinConstrainedProblem lcp, PBEproblem bound, PBEproblem unbound);
void LinConstrainedProblem_solve(LinConstrainedProblem lcp, Vector optimalCharges, Vector lagrangeMultipliers,
											Vector **approxCharges, int *numAdjustments, int* indexLooser, int* indexTighter);
void LinConstrainedProblem_getApproxSWMSoln(LinConstrainedProblem lcp, Vector newsoln, Vector oldsoln,
														  unsigned int orig, unsigned int newInd, real penalty);
void LinConstrainedProblem_operatorSave(LinConstrainedProblem lcp, char *filename);
void LinConstrainedProblem_preconditionerSave(LinConstrainedProblem lcp, char *filename);
// private interface: users should never have to screw with these
void LinConstrainedProblem_setupPreconditioner(LinConstrainedProblem lcp);
void LinConstrainedProblem_setupPreconditionerOrig(LinConstrainedProblem lcp);
void LinConstrainedProblem_setupRHS(LinConstrainedProblem lcp, Vector RHS);
void LinConstrainedProblem_GMRES(LinConstrainedProblem lcp, Vector soln, Vector RHS);
void LinConstrainedProblem_operatorMultiply(LinConstrainedProblem lcp, Vector Ax, Vector x);
void LinConstrainedProblem_preconditionerMultiply(LinConstrainedProblem lcp, Vector Px, Vector x); // new as of 2/28/07: block LU inversion
void LinConstrainedProblem_preconditionerMultiplyOrig(LinConstrainedProblem lcp, Vector Px, Vector x);
void LinConstrainedProblem_join(LinConstrainedProblem lcp, Vector dest, Vector q_temp, Vector lagrange_temp, Vector phi_b_temp, Vector phi_u_temp);
void LinConstrainedProblem_split(LinConstrainedProblem lcp, Vector q_temp, Vector lagrange_temp, Vector phi_b_temp, Vector phi_u_temp, Vector src);

#endif
