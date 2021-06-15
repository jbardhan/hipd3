#ifndef __BOXCONSTRAINED_H__
#define __BOXCONSTRAINED_H__

#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"
#include "Optimizer.h"
#include "Unconstrained.h"
#include "LinConstrained.h"

typedef struct _BoxConstrainedProblem {
   PBEproblem bound;
   PBEproblem unbound;

   Vector linearTerm;
   Vector Lqcur; // update at each iter
   
   unsigned int numconstraints;
   Matrix A_c;
   Vector b;
   
   Vector lowerBounds;
   Vector upperBounds;

   Matrix upperLeft;
   Preconditioner upperLeftP;

  Matrix penaltyMatrix;
  real penalty;

  Matrix leftSingularVectors;
  Matrix rightSingularVectors;
  Vector singularValues;
  
  real tolerance; 
  real looserTol;
  real tighterTol;

  Matrix Lhat;
} _BoxConstrainedProblem;

typedef struct _BoxConstrainedProblem* BoxConstrainedProblem;

// public interface
BoxConstrainedProblem BoxConstrainedProblem_allocate();
void BoxConstrainedProblem_free(BoxConstrainedProblem bcp);
void BoxConstrainedProblem_loadConstraints(BoxConstrainedProblem bcp, char *filename);
void BoxConstrainedProblem_setConstraints(BoxConstrainedProblem bcp,
                                          unsigned int numEqConstraints,
                                          Matrix A_c, Vector b,
                                          Vector lowerBounds, Vector upperBounds);
void BoxConstrainedProblem_setPBEproblems(BoxConstrainedProblem bcp,
                                          PBEproblem bound,
                                          PBEproblem unbound);
void BoxConstrainedProblem_setPenalty(BoxConstrainedProblem bcp, real penalty, Matrix penaltyMatrix,
												  Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues,
												  real tolerance, real looserTol, real tighterTol);

void BoxConstrainedProblem_solve(BoxConstrainedProblem bcp, Vector optimalCharges,
                                 Vector lagrangeMultipliers); // THERE IS NO SWM-EQUIV FOR BCP RIGHT NOW (3/29)
void BoxConstrainedProblem_operatorSave(BoxConstrainedProblem bcp, char *filename);
void BoxConstrainedProblem_preconditionerSave(BoxConstrainedProblem bcp, char *filename);

// private interface: user shouldn't have to screw with these
void BoxConstrainedProblem_setupUpperLeftBlock(BoxConstrainedProblem bcp, Vector q, Vector slack);
void BoxConstrainedProblem_calculateNewtonStep(BoxConstrainedProblem bcp, Vector step, Vector current, real sigma);
void BoxConstrainedProblem_setupRHS(BoxConstrainedProblem bcp, Vector RHS, Vector current, real sigma);
void BoxConstrainedProblem_GMRES(BoxConstrainedProblem bcp, Vector soln, Vector RHS); 
void BoxConstrainedProblem_operatorMultiply(BoxConstrainedProblem bcp, Vector Ax, Vector x);
void BoxConstrainedProblem_preconditionerMultiply(BoxConstrainedProblem bcp, Vector Px, Vector x);
void BoxConstrainedProblem_preconditionerMultiplyOrig(BoxConstrainedProblem bcp, Vector Px, Vector x);
real BoxConstrainedProblem_getStepScalingFactor(BoxConstrainedProblem bcp, Vector current, Vector step);
void BoxConstrainedProblem_update(BoxConstrainedProblem bcp, Vector current, Vector step, real alpha);
real BoxConstrainedProblem_checkConvergence(BoxConstrainedProblem bcp, Vector q, Vector slack);
void BoxConstrainedProblem_join(BoxConstrainedProblem bcp, Vector dest,
                                Vector q, Vector lagrange, Vector slack,
                                Vector phi_b, Vector phi_u);
void BoxConstrainedProblem_split(BoxConstrainedProblem bcp, 
                                Vector q, Vector lagrange, Vector slack,
                                Vector phi_b, Vector phi_u,
                                Vector src);
void BoxConstrainedProblem_computeDeltaReactionPotential(BoxConstrainedProblem bcp,
                                                         Vector Lq, Vector phi_b, Vector phi_u, Vector q);
#endif
