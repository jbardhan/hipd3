#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "FFTSVDpbeAPI.h"
#include "PBEproblem.h"

#define MAXITERATIONS 100
extern real tol;
extern int saveGMRES;
typedef enum { ABSOLUTE_TOL, RELATIVE_LHAT_TOL, RELATIVE_RAYLEIGH_TOL} penaltyToleranceType;
void Optimizer_computeHessian(PBEproblem bound, PBEproblem unbound, Matrix L, unsigned int startcol, unsigned endcol);
void Optimizer_computePreconditionerHessian(PBEproblem bound, PBEproblem unbound, Matrix Lhat);
void Optimizer_computeLinearTerm(PBEproblem bound, PBEproblem unbound, Vector c);
void Optimizer_zeroChainCharges(PBEproblem problem,  Vector globalCharges);
void Optimizer_getCoulombicInteraction(PBEproblem problem, Vector interaction);
void Optimizer_computePenaltyMatrix(PBEproblem bound, PBEproblem unbound, unsigned int problemsize, Matrix Lhat,
												penaltyToleranceType penaltyType, real tolerance, real penaltyScale, Matrix penalty,
												Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues, real *maxOpt);
void Optimizer_multiplyByL(PBEproblem bound, PBEproblem unbound, Vector Lx, Vector x);
void Optimizer_updateChargeDistributions(PBEproblem bound, PBEproblem unbound);
#endif
