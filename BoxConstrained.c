#include "BoxConstrained.h"

BoxConstrainedProblem BoxConstrainedProblem_allocate() {
   BoxConstrainedProblem problem;
   problem = (BoxConstrainedProblem)calloc(1, sizeof(_BoxConstrainedProblem));
   problem->numconstraints = 0;
   problem->bound = NULL;
   problem->unbound = NULL;
   problem->linearTerm = NULL;
   problem->A_c = NULL;
   problem->b = NULL;
   problem->lowerBounds = NULL;
   problem->upperBounds = NULL;
   problem->upperLeft = NULL;
   problem->upperLeftP = NULL;
	problem->Lhat = NULL;
	problem->penaltyMatrix = NULL;
   return problem;
}

void BoxConstrainedProblem_free(BoxConstrainedProblem bcp) {
   Vector_free(bcp->linearTerm);
   Matrix_free(bcp->A_c);
   Vector_free(bcp->b);
   Vector_free(bcp->lowerBounds);
   Vector_free(bcp->upperBounds);
   Matrix_free(bcp->upperLeft);
	Matrix_free(bcp->penaltyMatrix);
	Matrix_free(bcp->Lhat);
   Preconditioner_free(bcp->upperLeftP);
   free(bcp);
}

void BoxConstrainedProblem_setPenalty(BoxConstrainedProblem bcp, real penalty, Matrix penaltyMatrix,
												  Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues,
												  real tolerance, real looserTol, real tighterTol) {
  if (bcp->penaltyMatrix != NULL) {
	 Matrix_free(bcp->penaltyMatrix);
  }
  bcp->penalty = penalty;
  bcp->penaltyMatrix = Matrix_allocate(bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
  Matrix_copy(bcp->penaltyMatrix, penaltyMatrix, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
  bcp->tolerance = tolerance;
  bcp->looserTol = looserTol;
  bcp->tighterTol = tighterTol;
  if (bcp->rightSingularVectors != NULL) {
	 Matrix_free(bcp->rightSingularVectors);
  }
  bcp->rightSingularVectors = Matrix_allocate(bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
  Matrix_copy(bcp->rightSingularVectors, rightSingularVectors,
				  bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
  if (bcp->singularValues != NULL) {
	 Vector_free(bcp->singularValues);
  }
  bcp->singularValues = Vector_allocate(bcp->bound->numvariablecharges);
  Vector_copy(bcp->singularValues, singularValues, bcp->bound->numvariablecharges);
}

void BoxConstrainedProblem_loadConstraints(BoxConstrainedProblem bcp,
                                           char *filename) {
}

void BoxConstrainedProblem_setConstraints(BoxConstrainedProblem bcp,
                                          unsigned int numconstraints,
                                          Matrix A_c, Vector b,
                                          Vector lowerBounds, Vector upperBounds) {
   if (bcp->bound == NULL) {
      printf("BoxConstrainedProblem_setConstraints: set PBE problems before constraints!\n");
      exit(-1);
   }
   bcp->numconstraints = numconstraints;
	if (bcp->A_c != NULL) {
	  Matrix_free(bcp->A_c);
	}
   bcp->A_c = Matrix_allocate(numconstraints, bcp->bound->numvariablecharges);
   Matrix_copy(bcp->A_c, A_c, numconstraints, bcp->bound->numvariablecharges);
	if (bcp->b != NULL) {
	  Vector_free(bcp->b);
	}
	bcp->b = Vector_allocate(numconstraints);
   Vector_copy(bcp->b, b, numconstraints);
	if (bcp->lowerBounds != NULL) {
	  Vector_free(bcp->lowerBounds);
	  Vector_free(bcp->upperBounds);
	}
   bcp->lowerBounds = Vector_allocate(bcp->bound->numvariablecharges);
   bcp->upperBounds = Vector_allocate(bcp->bound->numvariablecharges);
   Vector_copy(bcp->lowerBounds, lowerBounds, bcp->bound->numvariablecharges);
   Vector_copy(bcp->upperBounds, upperBounds, bcp->bound->numvariablecharges);
}

void BoxConstrainedProblem_setPBEproblems(BoxConstrainedProblem bcp,
                                          PBEproblem bound, PBEproblem unbound) {
   bcp->bound = bound;
   bcp->unbound = unbound;
}

real BoxConstrainedProblem_checkConvergence(BoxConstrainedProblem bcp, Vector q, Vector slack) {
   return Vector_dot(q, slack, 2 * bcp->bound->numvariablecharges);   
}

real BoxConstrainedProblem_getStepScalingFactor(BoxConstrainedProblem bcp, Vector current, Vector step) {
   real alpha = 1.0;  
   unsigned int i;
   Vector cur_q, cur_lagrange, cur_slack, cur_phi_b, cur_phi_u;
   Vector step_q, step_lagrange, step_slack, step_phi_b, step_phi_u;
   cur_q = Vector_allocate(2 * bcp->bound->numvariablecharges);
   cur_lagrange = Vector_allocate(bcp->numconstraints + bcp->bound->numvariablecharges);
   cur_slack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   cur_phi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   cur_phi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   step_q = Vector_allocate(2 * bcp->bound->numvariablecharges);
   step_lagrange = Vector_allocate(bcp->numconstraints + bcp->bound->numvariablecharges);
   step_slack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   step_phi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   step_phi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

   BoxConstrainedProblem_split(bcp, cur_q, cur_lagrange, cur_slack, cur_phi_b, cur_phi_u, current);
   BoxConstrainedProblem_split(bcp, step_q, step_lagrange, step_slack, step_phi_b, step_phi_u, step);
   
   for (i = 0; i < 2 * bcp->bound->numvariablecharges; i++) {
      if ( ( step_q[i] < 0 ) && ( -cur_q[i] / step_q[i] < .999 * alpha) )         
         alpha = - .999 * cur_q[i] / step_q[i];
      if ( ( step_slack[i] < 0 ) && ( -cur_slack[i] / step_slack[i] < .999 * alpha) )
         alpha = - .999 * cur_slack[i] / step_slack[i];
   }

   Vector_free(cur_q);
   Vector_free(cur_lagrange);
   Vector_free(cur_slack);
   Vector_free(cur_phi_b);
   Vector_free(cur_phi_u);
   Vector_free(step_q);
   Vector_free(step_lagrange);
   Vector_free(step_slack);
   Vector_free(step_phi_b);
   Vector_free(step_phi_u);
   return alpha;
}

void BoxConstrainedProblem_solve(BoxConstrainedProblem bcp, Vector optimalCharges,
                                 Vector lagrangeMultipliers) {
   unsigned int iteration = 0;
   unsigned int maxIter = 50;
   Vector q, lagrange, slack, phi_b, phi_u;
   Vector qcur;
   Vector current, step;
   real alpha;
   FILE *OUT;
   unsigned int i;
   real convergenceTol = 1e-5 * bcp->bound->numvariablecharges;
   real convergenceCriterion = 1e8 * convergenceTol;
   real sigma = 0.4;
   qcur = Vector_allocate(bcp->bound->numvariablecharges);
   q = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange = Vector_allocate( bcp->bound->numvariablecharges + bcp->numconstraints);
   slack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   current = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints
                              + (bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables));
   step = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints
                              + (bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables));
   
   // initialize
   // step 0. set variables
   for (i = 0; i < 2 * bcp->bound->numvariablecharges; i++) {
      q[i] = 0.85;
      slack[i] = 1.0/ q[i];
   }
   Vector_zero(lagrange, bcp->bound->numvariablecharges + bcp->numconstraints);
   Vector_zero(phi_b, bcp->bound->numtotalsurfacevariables);
   Vector_zero(phi_u, bcp->unbound->numtotalsurfacevariables);
   Vector_copy(qcur, q, bcp->bound->numvariablecharges);
   Vector_addvector(qcur, bcp->lowerBounds, bcp->bound->numvariablecharges);

   // step 1. compute Lm
   bcp->Lqcur = Vector_allocate(bcp->bound->numvariablecharges);
   BoxConstrainedProblem_computeDeltaReactionPotential(bcp, bcp->Lqcur, phi_b, phi_u, qcur);
	// step 1a.) add in penalty term penalty*qcur;
	if (bcp->penaltyMatrix != NULL) {
	  Vector tmpPenaltyQ = Vector_allocate(bcp->bound->numvariablecharges);
	  Matrix_multiplyvector(tmpPenaltyQ, bcp->penaltyMatrix, qcur, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	  Vector_addvector(bcp->Lqcur, tmpPenaltyQ, bcp->bound->numvariablecharges);
	  Vector_free(tmpPenaltyQ);
	}
   // step 2.  compute linear term
   bcp->linearTerm = Vector_allocate(bcp->bound->numvariablecharges);
   Optimizer_computeLinearTerm(bcp->bound, bcp->unbound, bcp->linearTerm);
/*    OUT = fopen("linear","w"); */
/*    for (i = 0; i < bcp->bound->numvariablecharges; i++) */
/*       fprintf(OUT, "%f\n", bcp->linearTerm[i]); */
/*    fclose(OUT); */
   // end initialization

	unsigned int total_GMRES_iter = 0;
   BoxConstrainedProblem_join(bcp, current, q, lagrange, slack, phi_b, phi_u);
   while ( (convergenceCriterion > convergenceTol) && ( iteration < maxIter) ) {
      BoxConstrainedProblem_setupUpperLeftBlock(bcp, q, slack);
      BoxConstrainedProblem_calculateNewtonStep(bcp, step, current, sigma);
		total_GMRES_iter += num_GMRES_iter;
      alpha = BoxConstrainedProblem_getStepScalingFactor(bcp, current, step);
		printf("total GMRES so far = %d\n", total_GMRES_iter);
      BoxConstrainedProblem_update(bcp, current, step, alpha);
/* 		printf("total GMRES so far = %d\n", total_GMRES_iter); */
      BoxConstrainedProblem_split(bcp, q, lagrange, slack, phi_b, phi_u, current);
/* 		printf("total GMRES so far = %d\n", total_GMRES_iter); */
      convergenceCriterion = BoxConstrainedProblem_checkConvergence(bcp, q, slack);
      printf("Iteration %d:  alpha = %f   slackness violation = %f\n", iteration, alpha, convergenceCriterion);
      sigma = 0.4;
      iteration = iteration + 1;
   }

	num_GMRES_iter = total_GMRES_iter;
   if (iteration == maxIter) {
      printf("BoxConstrainedProblem_solve: failed to converge in %d Newton steps.\n", maxIter);
   } else {
      printf("BoxConstrainedProblem_solve: succeeded after %d Newton steps.\n", iteration);
   }

   Vector_copypiece(optimalCharges, 0, q, 0, bcp->bound->numvariablecharges);
   Vector_addvector(optimalCharges, bcp->lowerBounds, bcp->bound->numvariablecharges);
   Vector_copy(lagrangeMultipliers, lagrange, bcp->numconstraints);

/*    OUT = fopen("final", "w"); */
/*    for (i = 0; i < 5 * bcp->bound->numvariablecharges + bcp->numconstraints + bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables; i++) */
/*       fprintf(OUT, "%f\n", current[i]); */
/*    fclose(OUT); */

   printf("BoxConstrainedProblem_solve: pointwise slackness report: \n");
   for (i = 0; i < 2 * bcp->bound->numvariablecharges; i++) {
      printf("%f\n", q[i] * slack[i]);
   }
   Vector_free(current);
   Vector_free(step);
   Vector_free(q);
   Vector_free(lagrange);
   Vector_free(slack);
   Vector_free(phi_b);
   Vector_free(phi_u);
   Vector_free(qcur);
}

void BoxConstrainedProblem_update(BoxConstrainedProblem bcp, Vector current, Vector step, real alpha) {
   Vector dq, dlagrange, dslack, dphi_b, dphi_u, dLq;
   dq = Vector_allocate(2 * bcp->bound->numvariablecharges);
   dlagrange = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   dslack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   dphi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   dphi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   dLq = Vector_allocate(bcp->bound->numvariablecharges);
   BoxConstrainedProblem_split(bcp, dq, dlagrange, dslack, dphi_b, dphi_u, step);
   Vector_addscaledvector(current, alpha, step, 5 * bcp->bound->numvariablecharges + bcp->numconstraints
                          + (bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables));
   
   Vector_scale(dphi_b, alpha, bcp->bound->numtotalsurfacevariables);
   Vector_scale(dphi_u, alpha, bcp->unbound->numtotalsurfacevariables);
   PBEproblem_applyA3(bcp->bound, dLq, dphi_b);
   Vector_addvector(bcp->Lqcur, dLq, bcp->bound->numvariablecharges);
   PBEproblem_applyA3(bcp->unbound, dLq, dphi_u);
   Vector_addscaledvector(bcp->Lqcur, -1.0, dLq, bcp->bound->numvariablecharges);
	if (bcp->penaltyMatrix != NULL) {
	  Vector_zero(dLq, bcp->bound->numvariablecharges);
	  Matrix_multiplyvector(dLq, bcp->penaltyMatrix, dq, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	  Vector_addvector(bcp->Lqcur, dLq, bcp->bound->numvariablecharges);
	}
	
   Vector_free(dLq);
   Vector_free(dq);
   Vector_free(dlagrange);
   Vector_free(dslack);
   Vector_free(dphi_b);
   Vector_free(dphi_u);
}

void BoxConstrainedProblem_operatorSave(BoxConstrainedProblem bcp, char *filename) {
   unsigned int i, j;
   unsigned int matrix_size = 2 * bcp->bound->numvariablecharges + bcp->numconstraints 
      + bcp->bound->numvariablecharges + 2 * bcp->bound->numvariablecharges + bcp->bound->numtotalsurfacevariables
      + bcp->unbound->numtotalsurfacevariables;
   FILE *OUT = NULL;
   Vector x, Ax;
   x = Vector_allocate(matrix_size);
   Ax = Vector_allocate(matrix_size);
   OUT = fopen(filename, "w");
   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;
      BoxConstrainedProblem_operatorMultiply(bcp, Ax, x);
      for (j = 0; j < matrix_size; j++) 
         fprintf(OUT, "%f  ", Ax[j]);
      fprintf(OUT, "\n");
   }
   fclose(OUT);
   Vector_free(x);
   Vector_free(Ax);
}

void BoxConstrainedProblem_preconditionerSave(BoxConstrainedProblem bcp, char *filename) {
   unsigned int i, j;
   unsigned int matrix_size = 2 * bcp->bound->numvariablecharges + bcp->numconstraints + bcp->bound->numvariablecharges
      + 2 * bcp->bound->numvariablecharges + bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables;
   FILE *OUT = NULL;
   Vector x, Px;
   x = Vector_allocate(matrix_size);
   Px = Vector_allocate(matrix_size);
   OUT = fopen(filename, "w");
   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;
      BoxConstrainedProblem_preconditionerMultiply(bcp, Px, x);
      for (j = 0; j < matrix_size; j++) 
         fprintf(OUT, "%f  ", Px[j]);
      fprintf(OUT, "\n");
   }
   fclose(OUT);
   Vector_free(x);
   Vector_free(Px);
}

void BoxConstrainedProblem_setupUpperLeftBlock(BoxConstrainedProblem bcp,
                                               Vector q, Vector slack) {
   unsigned int block_size = 5 * bcp->bound->numvariablecharges + bcp->numconstraints;
   unsigned int i, j, index;
   if (bcp->upperLeft != NULL) { // will skip this only once (on the first Newton iteration)
      Matrix_free(bcp->upperLeft);
      Preconditioner_free(bcp->upperLeftP);
   }
   bcp->upperLeft = Matrix_allocate(block_size, block_size);
   bcp->upperLeftP = Preconditioner_allocate(block_size, block_size);

	if (bcp->Lhat == NULL) {
	  bcp->Lhat = Matrix_allocate(bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	  Optimizer_computePreconditionerHessian(bcp->bound, bcp->unbound, bcp->Lhat);
	}
	Matrix LhatPlusPenalty = Matrix_allocate(bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	Matrix_copy(LhatPlusPenalty, bcp->Lhat, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	
	// fill in penalty box if it exists
	if (bcp->penaltyMatrix != NULL) {
	  // moved the below up into this if()
	  Matrix_copypiece(bcp->upperLeft, 0, 0, bcp->penaltyMatrix, 0, 0, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	  Matrix_add(LhatPlusPenalty, LhatPlusPenalty, bcp->penaltyMatrix, bcp->bound->numvariablecharges, bcp->bound->numvariablecharges);
	}
	
	for (i = 0; i < bcp->bound->numvariablecharges; i++) {
	  for (j = 0; j < bcp->bound->numvariablecharges; j++) {
		 Preconditioner_set(bcp->upperLeftP, i, j, LhatPlusPenalty[i][j]);
	  }
	}
	Matrix_free(LhatPlusPenalty);
	
   // fill in Ac blocks
   for (i = 0; i < bcp->numconstraints; i++) {
      for (j = 0; j < bcp->bound->numvariablecharges; j++) {
         bcp->upperLeft[i+ 2 * bcp->bound->numvariablecharges][j] =
            bcp->A_c[i][j];
         bcp->upperLeft[j][2 * bcp->bound->numvariablecharges + i] =
            - bcp->A_c[i][j]; // NOTE THE MINUS SIGN!!
         Preconditioner_set(bcp->upperLeftP, i + 2 * bcp->bound->numvariablecharges, j, bcp->A_c[i][j]);
         Preconditioner_set(bcp->upperLeftP, j, i + 2 * bcp->bound->numvariablecharges, -bcp->A_c[i][j]);
      }
   }
   // fill in double primal I, I blocks
   for (i = 0; i < bcp->bound->numvariablecharges; i++) {
      bcp->upperLeft[i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges][i] = 1;
      bcp->upperLeft[i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges][i + bcp->bound->numvariablecharges] = 1;
      bcp->upperLeft[i][i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges] = -1;
      bcp->upperLeft[i + bcp->bound->numvariablecharges][i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges] = -1;
      Preconditioner_set(bcp->upperLeftP, i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges, i, 1.0);
      Preconditioner_set(bcp->upperLeftP, i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges, i + bcp->bound->numvariablecharges, 1.0);
      Preconditioner_set(bcp->upperLeftP, i, i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges, -1.0);
      Preconditioner_set(bcp->upperLeftP, i + bcp->bound->numvariablecharges, i + bcp->numconstraints + 2 * bcp->bound->numvariablecharges, -1.0);
   }

   // three square, diagonal blocks
   for (i = 0; i < 2 * bcp->bound->numvariablecharges; i++) {
      index = 3 * bcp->bound->numvariablecharges + bcp->numconstraints + i;
      // fill in dual slack I block
      bcp->upperLeft[i][index] = -1.0;
      Preconditioner_set(bcp->upperLeftP, i, index, -1.0);
      // fill in slack diagonal block 
      bcp->upperLeft[index][i] = slack[i];
      Preconditioner_set(bcp->upperLeftP, index, i, slack[i]);
      // fill in primal diagonal block
      bcp->upperLeft[index][index] = q[i];
      Preconditioner_set(bcp->upperLeftP, index, index, q[i]);
   }

   Preconditioner_factor(bcp->upperLeftP);
}

void BoxConstrainedProblem_calculateNewtonStep(BoxConstrainedProblem bcp,
                                               Vector step, Vector current, real sigma) {
   FILE *OUT;
   unsigned int i;
	static unsigned int iter = 0;
   unsigned int matrix_size = 5 * bcp->bound->numvariablecharges + bcp->numconstraints 
      + (bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables);
   Vector RHS = Vector_allocate(matrix_size);

   BoxConstrainedProblem_setupRHS(bcp, RHS, current, sigma);
   BoxConstrainedProblem_GMRES(bcp, step, RHS);

/* 	char filename[100]; */
/* 	sprintf(filename, "Abox_%d", iter); */
/*    BoxConstrainedProblem_operatorSave(bcp, filename); */
/* 	sprintf(filename, "Pbox_%d", iter); */
/*    BoxConstrainedProblem_preconditionerSave(bcp, filename); */

/* 	Matrix_writefile("uleft", bcp->upperLeft, 5 * bcp->bound->numvariablecharges + bcp->numconstraints , 5 * bcp->bound->numvariablecharges + bcp->numconstraints  ); */
/* 	OUT = fopen("rhs","w"); */
/*    for (i = 0; i < matrix_size; i++) */
/*       fprintf(OUT, "%f\n", RHS[i]); */
/*    fclose(OUT); */
/*    OUT = fopen("current", "w"); */
/*    for (i = 0; i < matrix_size; i++) */
/*       fprintf(OUT, "%f\n", current[i]); */
/*    fclose(OUT); */
/*    OUT = fopen("step", "w"); */
/*    for (i = 0; i < matrix_size; i++) */
/*       fprintf(OUT, "%f\n", step[i]); */
/*    fclose(OUT); */
/* 	exit(-1); */
	iter = iter + 1;
	
   Vector_free(RHS);
}

void BoxConstrainedProblem_setupRHS(BoxConstrainedProblem bcp, Vector RHS, Vector current, real sigma) {
                                    
   unsigned int i;
   real gamma;
   Vector q, lagrange, slack, phi_b, phi_u;
   Vector rhs_q, rhs_lagrange, rhs_slack, rhs_phi_b, rhs_phi_u;
   Vector rhs_q_work, rhs_lagrange_work, rhs_slack_work, rhs_phi_b_work, rhs_phi_u_work;
   q = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange = Vector_allocate(bcp->numconstraints + bcp->bound->numvariablecharges);
   slack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   rhs_q = Vector_allocate(2 * bcp->bound->numvariablecharges);
   rhs_lagrange = Vector_allocate(bcp->numconstraints + bcp->bound->numvariablecharges);
   rhs_slack = Vector_allocate(2 * bcp->bound->numvariablecharges);
   rhs_phi_b = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   rhs_phi_u = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   rhs_q_work = Vector_allocate(2 * bcp->bound->numvariablecharges);
   rhs_lagrange_work = Vector_allocate(bcp->numconstraints + bcp->bound->numvariablecharges);
   rhs_slack_work = Vector_allocate(2 * bcp->bound->numvariablecharges);
   rhs_phi_b_work = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   rhs_phi_u_work = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   
   BoxConstrainedProblem_split(bcp, q, lagrange, slack, phi_b, phi_u, current);
   
   Vector_addscaledvector(rhs_q, -1.0, bcp->linearTerm, bcp->bound->numvariablecharges);
   Vector_addscaledvector(rhs_q, -1.0, bcp->Lqcur, bcp->bound->numvariablecharges); // don't forget to update Lqcur each iter!!
   Vector_addscaledvector(rhs_q, 1.0, slack, 2 * bcp->bound->numvariablecharges); // this adds to both
   Matrix_multiplyvector_transpose(rhs_q_work, bcp->A_c, lagrange, bcp->numconstraints, bcp->bound->numvariablecharges);
   Vector_addvector(rhs_q, rhs_q_work, bcp->bound->numvariablecharges);
   Vector_addvector(rhs_q, &(lagrange[bcp->numconstraints]), bcp->bound->numvariablecharges);
   Vector_addvector(&(rhs_q[bcp->bound->numvariablecharges]), &(lagrange[bcp->numconstraints]), bcp->bound->numvariablecharges);

   Vector_addvector(rhs_lagrange, bcp->b, bcp->numconstraints);
   
   Matrix_multiplyvector(rhs_lagrange_work, bcp->A_c, bcp->lowerBounds, bcp->numconstraints, bcp->bound->numvariablecharges);
   Vector_addscaledvector(rhs_lagrange, -1.0, rhs_lagrange_work, bcp->numconstraints);
   Vector_zero(rhs_lagrange_work, bcp->numconstraints);
   Matrix_multiplyvector(rhs_lagrange_work, bcp->A_c, q, bcp->numconstraints, bcp->bound->numvariablecharges);
   Vector_addscaledvector(rhs_lagrange, -1.0, rhs_lagrange_work, bcp->numconstraints);
   // second half of rhs_lagrange
   Vector_addvector(&(rhs_lagrange[bcp->numconstraints]), bcp->upperBounds, bcp->bound->numvariablecharges);
   Vector_addscaledvector(&(rhs_lagrange[bcp->numconstraints]), -1.0, bcp->lowerBounds, bcp->bound->numvariablecharges);
   Vector_addscaledvector(&(rhs_lagrange[bcp->numconstraints]), -1.0, q, bcp->bound->numvariablecharges);
   Vector_addscaledvector(&(rhs_lagrange[bcp->numconstraints]), -1.0, &(q[bcp->bound->numvariablecharges]), bcp->bound->numvariablecharges);

   gamma = Vector_dot(q, slack, 2 * bcp->bound->numvariablecharges);
   gamma = gamma / (2 * bcp->bound->numvariablecharges);
   for (i = 0; i < 2 * bcp->bound->numvariablecharges; i++)
      rhs_slack[i] = -q[i] * slack[i] + sigma * gamma; // centering parameter shows up here
   
   PBEproblem_applyA1(bcp->bound, rhs_phi_b_work, bcp->lowerBounds);
   Vector_addscaledvector(rhs_phi_b, +1.0, rhs_phi_b_work, bcp->bound->numtotalsurfacevariables);
   Vector_zero(rhs_phi_b_work, bcp->bound->numtotalsurfacevariables);
   PBEproblem_applyA1(bcp->bound, rhs_phi_b_work, q);
   Vector_addscaledvector(rhs_phi_b, +1.0, rhs_phi_b_work, bcp->bound->numtotalsurfacevariables);
   Vector_zero(rhs_phi_b_work, bcp->bound->numtotalsurfacevariables);
   PBEproblem_applyA2(bcp->bound, rhs_phi_b_work, phi_b);
   Vector_addscaledvector(rhs_phi_b, -1.0, rhs_phi_b_work, bcp->bound->numtotalsurfacevariables);
   
   PBEproblem_applyA1(bcp->unbound, rhs_phi_u_work, q);
   Vector_addscaledvector(rhs_phi_u, +1.0, rhs_phi_u_work, bcp->unbound->numtotalsurfacevariables);
   Vector_zero(rhs_phi_u_work, bcp->unbound->numtotalsurfacevariables);
   PBEproblem_applyA1(bcp->unbound, rhs_phi_u_work, bcp->lowerBounds);
   Vector_addscaledvector(rhs_phi_u, +1.0, rhs_phi_u_work, bcp->unbound->numtotalsurfacevariables);
   Vector_zero(rhs_phi_u_work, bcp->unbound->numtotalsurfacevariables);
   PBEproblem_applyA2(bcp->unbound, rhs_phi_u_work, phi_u);
   Vector_addscaledvector(rhs_phi_u, -1.0, rhs_phi_u_work, bcp->unbound->numtotalsurfacevariables);

   BoxConstrainedProblem_join(bcp, RHS, rhs_q, rhs_lagrange, rhs_slack, rhs_phi_b, rhs_phi_u);

   Vector_free(rhs_q_work);
   Vector_free(rhs_lagrange_work);
   Vector_free(rhs_slack_work);
   Vector_free(rhs_phi_b_work);
   Vector_free(rhs_phi_u_work);
   Vector_free(rhs_q);
   Vector_free(rhs_lagrange);
   Vector_free(rhs_slack);
   Vector_free(rhs_phi_b);
   Vector_free(rhs_phi_u);
   Vector_free(q);
   Vector_free(lagrange);
   Vector_free(slack);
   Vector_free(phi_b);
   Vector_free(phi_u);
}

void BoxConstrainedProblem_GMRES(BoxConstrainedProblem bcp, Vector sol, Vector rhs) {
   unsigned int size = 5 * bcp->bound->numvariablecharges + bcp->numconstraints
      + (bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables);
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   BoxConstrainedProblem_preconditionerMultiply(bcp, r, rhs);

   normr = Vector_norm(r, size);

   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);

   for (i = 0; i < MAXITERATIONS; i++) {
      BoxConstrainedProblem_operatorMultiply(bcp, P, bv[i]);
      BoxConstrainedProblem_preconditionerMultiply(bcp, P, P);
      
      for (j = 0; j <= i; j++) {
         H[i][j] = Vector_dot(P, bv[j], size);
         Vector_subtractscaledvector(P, H[i][j], bv[j], size);
      }
      
      H[i][i+1] = Vector_norm(P, size);
      bv[i+1] = Vector_allocate(size);
      Vector_copy(bv[i+1], P, size);
      Vector_scale(bv[i+1], 1.0 / H[i][i+1], size);
      
      for (k = 0; k < i; k++)
         givensrotate(c[k], s[k], &H[i][k], &H[i][k+1]);
      givens(H[i][i], H[i][i+1], &c[i], &s[i]);
      givensrotate(c[i], s[i], &H[i][i], &H[i][i+1]);
      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %2.8f\n", i+1, residual);
#endif
      if (residual < tol)
         break;
   }
   for (k = 0; k <= i; k++)
      y[k] = g[k];

	num_GMRES_iter = i;
      
   for (k = i; k >= 0; k--) {
      y[k] /= H[k][k];
      for (j = k-1; j >= 0; j--)
         y[j] -= H[k][j] * y[k];
   }

   for (j = 0; j <= i; j++)
      Vector_addscaledvector(x, y[j], bv[j], size);

   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(r);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);

   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);   
}

void BoxConstrainedProblem_operatorMultiply(BoxConstrainedProblem bcp, Vector Ax, Vector x) {
   unsigned int block_size = 5 * bcp->bound->numvariablecharges + bcp->numconstraints;
   Vector chargeAndLagrangeAndSlack, AchargeAndLagrangeAndSlack;
   Vector q_temp, lagrange_temp, slack_temp, phi_b_temp, phi_u_temp;
   Vector q_work, lagrange_work, slack_work, phi_b_work, phi_u_work;
   Vector Aq_temp, Alagrange_temp, Aslack_temp, Aphi_b_temp, Aphi_u_temp;

   chargeAndLagrangeAndSlack = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   AchargeAndLagrangeAndSlack = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   q_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   slack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

   q_work = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange_work = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   slack_work = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b_work = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u_work = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

   Aq_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Alagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   Aslack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Aphi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   Aphi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   
   BoxConstrainedProblem_split(bcp, q_temp, lagrange_temp, slack_temp,
                               phi_b_temp, phi_u_temp, x);

   // BEGIN MULTIPLY STUFF

   //  ROWS 1, 2, 3: upperLeft
   Vector_copypiece(chargeAndLagrangeAndSlack, 0, q_temp, 0,
                    2 * bcp->bound->numvariablecharges);
   Vector_copypiece(chargeAndLagrangeAndSlack, 2 * bcp->bound->numvariablecharges,
                    lagrange_temp, 0,
                    bcp->bound->numvariablecharges + bcp->numconstraints);
   Vector_copypiece(chargeAndLagrangeAndSlack, 3 * bcp->bound->numvariablecharges + bcp->numconstraints,
                    slack_temp, 0,
                    2 * bcp->bound->numvariablecharges); // penalty is already in bcp->upperLeft!
   Matrix_multiplyvector(AchargeAndLagrangeAndSlack, bcp->upperLeft, chargeAndLagrangeAndSlack, block_size, block_size);
   Vector_copypiece(Aq_temp, 0, AchargeAndLagrangeAndSlack, 0,
                    2 * bcp->bound->numvariablecharges);
   Vector_copypiece(Alagrange_temp, 0, AchargeAndLagrangeAndSlack,
                    2 * bcp->bound->numvariablecharges,
                    bcp->bound->numvariablecharges + bcp->numconstraints);
   Vector_copypiece(Aslack_temp, 0, AchargeAndLagrangeAndSlack,
                    3 * bcp->bound->numvariablecharges + bcp->numconstraints,
                    2 * bcp->bound->numvariablecharges);

   //  ROW 1: A3b, -A3u.  notice that i'm cheating by adding only the first numvarcharges
   //  (recall we have to expand the box constraints so we have 2 n_c primal vars)
   PBEproblem_applyA3(bcp->bound, q_work, phi_b_temp);
   Vector_addvector(Aq_temp, q_work, bcp->bound->numvariablecharges);
   PBEproblem_applyA3(bcp->unbound, q_work, phi_u_temp);
   Vector_addscaledvector(Aq_temp, -1.0, q_work, bcp->bound->numvariablecharges); // subtracting Lb-Lu!
   
   //  ROW 4: -A1b, A2b
   PBEproblem_applyA1(bcp->bound, phi_b_work, q_temp);
   Vector_addscaledvector(Aphi_b_temp, -1.0, phi_b_work, bcp->bound->numtotalsurfacevariables); 
   PBEproblem_applyA2(bcp->bound, phi_b_work, phi_b_temp);
   Vector_addvector(Aphi_b_temp, phi_b_work, bcp->bound->numtotalsurfacevariables);

   //  ROW 5: -A1u, A2u
   PBEproblem_applyA1(bcp->unbound, phi_u_work, q_temp);
   Vector_addscaledvector(Aphi_u_temp, -1.0, phi_u_work, bcp->unbound->numtotalsurfacevariables); 
   PBEproblem_applyA2(bcp->unbound, phi_u_work, phi_u_temp);
   Vector_addvector(Aphi_u_temp, phi_u_work, bcp->unbound->numtotalsurfacevariables);
   
   // END MULTIPLY STUFF
   
   BoxConstrainedProblem_join(bcp, Ax, Aq_temp, Alagrange_temp, Aslack_temp,
                              Aphi_b_temp, Aphi_u_temp);

   Vector_free(chargeAndLagrangeAndSlack);
   Vector_free(AchargeAndLagrangeAndSlack);
   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(slack_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(q_work);
   Vector_free(lagrange_work);
   Vector_free(slack_work);
   Vector_free(phi_b_work);
   Vector_free(phi_u_work);
   Vector_free(Aq_temp);
   Vector_free(Alagrange_temp);
   Vector_free(Aslack_temp);
   Vector_free(Aphi_b_temp);
   Vector_free(Aphi_u_temp);

}

void BoxConstrainedProblem_preconditionerMultiply(BoxConstrainedProblem bcp, Vector Px, Vector x) {
   Vector chargeAndLagrange_temp, PchargeAndLagrange;
   Vector q_temp, lagrange_temp, slack_temp, phi_b_temp, phi_u_temp;
   Vector Pq_temp, Plagrange_temp, Pslack_temp, Pphi_b_temp, Pphi_u_temp;

	// identity preconditioner should be commented out just about all the time!!
/* 	Vector_copy(Px, x, 5 * bcp->bound->numvariablecharges + bcp->numconstraints + bcp->bound->numtotalsurfacevariables + bcp->unbound->numtotalsurfacevariables); */
/* 	return; */
	// end identity preconditioner
	
   chargeAndLagrange_temp = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   PchargeAndLagrange = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   q_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   slack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

   Pq_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Plagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   Pslack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Pphi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   Pphi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   
   BoxConstrainedProblem_split(bcp, q_temp, lagrange_temp, slack_temp,
                               phi_b_temp, phi_u_temp, x);


	unsigned int sizePD = 5 * bcp->bound->numvariablecharges + bcp->numconstraints;
	Vector P1chargeAndLagrange_temp, P1phi_b_temp, P1phi_u_temp;
	P1chargeAndLagrange_temp = Vector_allocate(sizePD);
	P1phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
	P1phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

	Vector P2P1chargeAndLagrange_temp, P2P1phi_b_temp, P2P1phi_u_temp;
	P2P1chargeAndLagrange_temp = Vector_allocate(sizePD);
	P2P1phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
	P2P1phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

	Vector P3P2P1chargeAndLagrange_temp, P3P2P1phi_b_temp, P3P2P1phi_u_temp;
	P3P2P1chargeAndLagrange_temp = Vector_allocate(sizePD);
	P3P2P1phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
	P3P2P1phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

	// BEGIN MULTIPLY STUFF

	// P1: [I 0 0; 0 Pb 0; 0 0 Pu]
	Vector_copypiece(P1chargeAndLagrange_temp, 0, q_temp, 0, 2 * bcp->bound->numvariablecharges);
	Vector_copypiece(P1chargeAndLagrange_temp, 2*bcp->bound->numvariablecharges, lagrange_temp, 0, bcp->numconstraints + bcp->bound->numvariablecharges);
	Vector_copypiece(P1chargeAndLagrange_temp, 3*bcp->bound->numvariablecharges+bcp->numconstraints,
						  slack_temp, 0, 2 * bcp->bound->numvariablecharges);

	PBEproblem_applyPreconditioner(bcp->bound, P1phi_b_temp, phi_b_temp);
	PBEproblem_applyPreconditioner(bcp->unbound, P1phi_u_temp, phi_u_temp);

	// P2: [I -A3b -A3u; 0 I 0; 0 0 I];
	Vector A3P1phi_b, A3P1phi_u;
	A3P1phi_b = Vector_allocate(bcp->bound->numvariablecharges);
	A3P1phi_u = Vector_allocate(bcp->bound->numvariablecharges);
	PBEproblem_applyA3(bcp->bound, A3P1phi_b, P1phi_b_temp); 
	PBEproblem_applyA3(bcp->unbound, A3P1phi_u, P1phi_u_temp);
	Vector_scale(A3P1phi_b, -1., bcp->bound->numvariablecharges);
	Vector_addvector(A3P1phi_b, A3P1phi_u, bcp->bound->numvariablecharges);
	Vector_copy(P2P1chargeAndLagrange_temp, P1chargeAndLagrange_temp, sizePD);
	Vector_addvector(P2P1chargeAndLagrange_temp, A3P1phi_b, bcp->bound->numvariablecharges);
	Vector_copy(P2P1phi_b_temp, P1phi_b_temp, bcp->bound->numtotalsurfacevariables); 
	Vector_copy(P2P1phi_u_temp, P1phi_u_temp, bcp->unbound->numtotalsurfacevariables);
	Vector_free(A3P1phi_b); Vector_free(A3P1phi_u);
	
	// P3: [upperLeftP 0 0; 0 I 0; 0 0 I];
	Preconditioner_solve(P3P2P1chargeAndLagrange_temp, bcp->upperLeftP, P2P1chargeAndLagrange_temp);
	Vector_copy(P3P2P1phi_b_temp, P2P1phi_b_temp, bcp->bound->numtotalsurfacevariables);
	Vector_copy(P3P2P1phi_u_temp, P2P1phi_u_temp, bcp->unbound->numtotalsurfacevariables);
	
	// P4: [I 0 0; Pb*A1b I 0; Pu *A1u 0 I]
	Vector A1bq_temp, A1uq_temp;
	A1bq_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
	A1uq_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
	Vector_copy(PchargeAndLagrange, P3P2P1chargeAndLagrange_temp, sizePD);//
	Vector_copypiece(q_temp, 0, P3P2P1chargeAndLagrange_temp, 0, bcp->bound->numvariablecharges);//
	PBEproblem_applyA1(bcp->bound, A1bq_temp, q_temp);
	PBEproblem_applyA1(bcp->unbound, A1uq_temp, q_temp);
	PBEproblem_applyPreconditioner(bcp->bound, Pphi_b_temp, A1bq_temp);
	PBEproblem_applyPreconditioner(bcp->unbound, Pphi_u_temp, A1uq_temp);
	Vector_addvector(Pphi_b_temp, P3P2P1phi_b_temp, bcp->bound->numtotalsurfacevariables);
	Vector_addvector(Pphi_u_temp, P3P2P1phi_u_temp, bcp->unbound->numtotalsurfacevariables);
	Vector_free(A1bq_temp);
	Vector_free(A1uq_temp);
	
   // END MULTIPLY STUFF
   
	Vector_copypiece(Pq_temp, 0, PchargeAndLagrange, 0, 2 * bcp->bound->numvariablecharges);
	Vector_copypiece(Plagrange_temp, 0, PchargeAndLagrange, 2 * bcp->bound->numvariablecharges, bcp->numconstraints + bcp->bound->numvariablecharges);
	Vector_copypiece(Pslack_temp, 0, PchargeAndLagrange, 3 * bcp->bound->numvariablecharges + bcp->numconstraints, 2 * bcp->bound->numvariablecharges);
	
   BoxConstrainedProblem_join(bcp, Px, Pq_temp, Plagrange_temp, Pslack_temp,
                              Pphi_b_temp, Pphi_u_temp);

	Vector_free(P1chargeAndLagrange_temp); Vector_free(P1phi_b_temp); Vector_free(P1phi_u_temp);
	Vector_free(P2P1chargeAndLagrange_temp); Vector_free(P2P1phi_b_temp); Vector_free(P2P1phi_u_temp);
	Vector_free(P3P2P1chargeAndLagrange_temp); Vector_free(P3P2P1phi_b_temp); Vector_free(P3P2P1phi_u_temp);

   Vector_free(chargeAndLagrange_temp);
   Vector_free(PchargeAndLagrange);
   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(slack_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Pq_temp);
   Vector_free(Plagrange_temp);
   Vector_free(Pslack_temp);
   Vector_free(Pphi_b_temp);
   Vector_free(Pphi_u_temp);
}

void BoxConstrainedProblem_preconditionerMultiplyOrig(BoxConstrainedProblem bcp, Vector Px, Vector x) {
   Vector chargeAndLagrangeAndSlack, PchargeAndLagrangeAndSlack;
   Vector q_temp, lagrange_temp, slack_temp, phi_b_temp, phi_u_temp;
   Vector Pq_temp, Plagrange_temp, Pslack_temp, Pphi_b_temp, Pphi_u_temp;

   chargeAndLagrangeAndSlack = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   PchargeAndLagrangeAndSlack = Vector_allocate( 5 * bcp->bound->numvariablecharges + bcp->numconstraints);
   q_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   lagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   slack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   phi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);

   Pq_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Plagrange_temp = Vector_allocate(bcp->bound->numvariablecharges + bcp->numconstraints);
   Pslack_temp = Vector_allocate(2 * bcp->bound->numvariablecharges);
   Pphi_b_temp = Vector_allocate(bcp->bound->numtotalsurfacevariables);
   Pphi_u_temp = Vector_allocate(bcp->unbound->numtotalsurfacevariables);
   
   BoxConstrainedProblem_split(bcp, q_temp, lagrange_temp, slack_temp,
                               phi_b_temp, phi_u_temp, x);

   // BEGIN MULTIPLY STUFF
   // ROWS 1, 2, and 3: Preconditioner for upper left block
   Vector_copypiece(chargeAndLagrangeAndSlack, 0, q_temp, 0,
                    2 * bcp->bound->numvariablecharges);
   Vector_copypiece(chargeAndLagrangeAndSlack, 2 * bcp->bound->numvariablecharges,
                    lagrange_temp, 0,
                    bcp->bound->numvariablecharges + bcp->numconstraints);
   Vector_copypiece(chargeAndLagrangeAndSlack, 3 * bcp->bound->numvariablecharges + bcp->numconstraints,
                    slack_temp, 0,
                    2 * bcp->bound->numvariablecharges);
   Preconditioner_solve(PchargeAndLagrangeAndSlack, bcp->upperLeftP, chargeAndLagrangeAndSlack);
   Vector_copypiece(Pq_temp, 0, PchargeAndLagrangeAndSlack, 0,
                    2 * bcp->bound->numvariablecharges);
   Vector_copypiece(Plagrange_temp, 0, PchargeAndLagrangeAndSlack,
                    2 * bcp->bound->numvariablecharges,
                    bcp->bound->numvariablecharges + bcp->numconstraints);
   Vector_copypiece(Pslack_temp, 0, PchargeAndLagrangeAndSlack,
                    3 * bcp->bound->numvariablecharges + bcp->numconstraints,
                    2 * bcp->bound->numvariablecharges);

   // ROW 4: A2b
   Preconditioner_solve(Pphi_b_temp, bcp->bound->preconditioner, phi_b_temp);
   // ROW 5: A2u
   Preconditioner_solve(Pphi_u_temp, bcp->unbound->preconditioner, phi_u_temp);
   // END MULTIPLY STUFF
   
   BoxConstrainedProblem_join(bcp, Px, Pq_temp, Plagrange_temp, Pslack_temp,
                              Pphi_b_temp, Pphi_u_temp);

   Vector_free(chargeAndLagrangeAndSlack);
   Vector_free(PchargeAndLagrangeAndSlack);
   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(slack_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Pq_temp);
   Vector_free(Plagrange_temp);
   Vector_free(Pslack_temp);
   Vector_free(Pphi_b_temp);
   Vector_free(Pphi_u_temp);

}

void BoxConstrainedProblem_join(BoxConstrainedProblem bcp, Vector dest,
                                Vector q, Vector lagrange, Vector slack,
                                Vector phi_b, Vector phi_u) {
   unsigned int cur = 0;
   Vector_copypiece(dest, cur, q, 0, 2 * bcp->bound->numvariablecharges);
   cur += 2 * bcp->bound->numvariablecharges;
   Vector_copypiece(dest, cur, lagrange, 0, bcp->bound->numvariablecharges + bcp->numconstraints);
   cur += bcp->bound->numvariablecharges + bcp->numconstraints;
   Vector_copypiece(dest, cur, slack, 0, 2 * bcp->bound->numvariablecharges);
   cur += 2 * bcp->bound->numvariablecharges;
   Vector_copypiece(dest, cur, phi_b, 0, bcp->bound->numtotalsurfacevariables);
   cur += bcp->bound->numtotalsurfacevariables;
   Vector_copypiece(dest, cur, phi_u, 0, bcp->unbound->numtotalsurfacevariables);
   cur += bcp->unbound->numtotalsurfacevariables;
}

void BoxConstrainedProblem_split(BoxConstrainedProblem bcp, 
                                Vector q, Vector lagrange, Vector slack,
                                Vector phi_b, Vector phi_u,
                                 Vector src) {
   unsigned int cur = 0;
   Vector_copypiece(q, 0, src, cur, 2* bcp->bound->numvariablecharges);
   cur += 2 * bcp->bound->numvariablecharges;
   Vector_copypiece(lagrange, 0, src, cur, bcp->bound->numvariablecharges + bcp->numconstraints);
   cur += bcp->bound->numvariablecharges + bcp->numconstraints;
   Vector_copypiece(slack, 0, src, cur, 2 * bcp->bound->numvariablecharges);
   cur += 2 * bcp->bound->numvariablecharges;
   Vector_copypiece(phi_b, 0, src, cur, bcp->bound->numtotalsurfacevariables);
   cur += bcp->bound->numtotalsurfacevariables;
   Vector_copypiece(phi_u, 0, src, cur, bcp->unbound->numtotalsurfacevariables);
}

void BoxConstrainedProblem_computeDeltaReactionPotential(BoxConstrainedProblem bcp,
                                                         Vector Lq, Vector phi_b, Vector phi_u, Vector q)
{
   Vector unboundLm, boundLm;
   unboundLm = Vector_allocate(bcp->bound->numvariablecharges);
   boundLm = Vector_allocate(bcp->bound->numvariablecharges);
   PBEproblem_setVariableChargeVector(bcp->unbound, q);
   PBEproblem_solve(bcp->unbound);
   PBEproblem_getVariableReactionPotentials(bcp->unbound, unboundLm);
   Vector_copy(phi_u, bcp->unbound->Sol, bcp->unbound->numtotalsurfacevariables);
   
   PBEproblem_setVariableChargeVector(bcp->bound, q);
   PBEproblem_solve(bcp->bound);
   PBEproblem_getVariableReactionPotentials(bcp->bound, boundLm);
   Vector_copy(phi_b, bcp->bound->Sol, bcp->bound->numtotalsurfacevariables);
   
   Vector_copy(Lq, boundLm, bcp->bound->numvariablecharges);
   Vector_addscaledvector(Lq, -1.0, unboundLm, bcp->bound->numvariablecharges);
   Vector_free(unboundLm);
   Vector_free(boundLm);
}
