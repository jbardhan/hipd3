#include "LinConstrained.h"
LinConstrainedProblem LinConstrainedProblem_allocate() {
   LinConstrainedProblem lcp;

   lcp = (LinConstrainedProblem)calloc(1, sizeof(_LinConstrainedProblem));
   lcp->A_c = NULL;
   lcp->b = NULL;
	lcp->Lhat = NULL;
   lcp->LhatAndAcInv = NULL;
   lcp->linearTerm = NULL;
	lcp->penaltyMatrix = NULL;
   return lcp;
}

void LinConstrainedProblem_free(LinConstrainedProblem lcp) {
   Matrix_free(lcp->A_c);
   Vector_free(lcp->b);
	Matrix_free(lcp->penaltyMatrix);
	Matrix_free(lcp->Lhat);
   Matrix_free(lcp->LhatAndAcInv);
   Vector_free(lcp->linearTerm);
   free(lcp);
}

void LinConstrainedProblem_setConstraints(LinConstrainedProblem lcp, unsigned int numconstraints,
                                          Matrix A_c, Vector b) {
  if (lcp->A_c != NULL) {
	 Matrix_free(lcp->A_c);
  }
  if (lcp->b != NULL) {
	 Vector_free(lcp->b);
  }

  lcp->numconstraints = numconstraints;
  lcp->A_c = Matrix_allocate(numconstraints, lcp->bound->numvariablecharges);
  Matrix_copy(lcp->A_c, A_c, numconstraints, lcp->bound->numvariablecharges);

  lcp->b   = Vector_allocate(numconstraints);
  Vector_copy(lcp->b, b, numconstraints);

  if (lcp->LhatAndAcInv != NULL) { // so that we reinitialize it on next solve
	 Matrix_free(lcp->LhatAndAcInv);
	 lcp->LhatAndAcInv = NULL;
  }
}

void LinConstrainedProblem_setPenalty(LinConstrainedProblem lcp, real penalty, Matrix penaltyMatrix,
												  Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues,
												  real tolerance, real looserTol, real tighterTol) {
  if (lcp->penaltyMatrix != NULL) {
	 Matrix_free(lcp->penaltyMatrix);
  }
  lcp->penalty = penalty;
  lcp->penaltyMatrix = Matrix_allocate(lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
  Matrix_copy(lcp->penaltyMatrix, penaltyMatrix, lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);

  if (lcp->LhatAndAcInv != NULL) { // so that we reinitialize it on next solve
	 Matrix_free(lcp->LhatAndAcInv);
	 lcp->LhatAndAcInv = NULL;
  }
}


void LinConstrainedProblem_setPBEproblems(LinConstrainedProblem lcp, PBEproblem bound, PBEproblem unbound) {
   lcp->bound = bound;
   lcp->unbound = unbound;
}

void LinConstrainedProblem_setupPreconditioner(LinConstrainedProblem lcp) {
  Matrix Lhatlocal, LhatAndAc,Actranspose;
  
  if (lcp->Lhat == NULL) {
	 lcp->Lhat = Matrix_allocate(lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
	 Optimizer_computePreconditionerHessian(lcp->bound, lcp->unbound, lcp->Lhat);
  }
  if (lcp->LhatAndAcInv != NULL)
	 Matrix_free(lcp->LhatAndAcInv);
  
  Lhatlocal = Matrix_allocate(lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
  lcp->LhatAndAcInv = Matrix_allocate(lcp->bound->numvariablecharges + lcp->numconstraints,
												  lcp->bound->numvariablecharges + lcp->numconstraints);
  LhatAndAc         = Matrix_allocate(lcp->bound->numvariablecharges + lcp->numconstraints,
                                       lcp->bound->numvariablecharges + lcp->numconstraints);
  Actranspose = Matrix_allocate(lcp->numconstraints, lcp->bound->numvariablecharges);
   
  Matrix_copy(Actranspose, lcp->A_c, lcp->numconstraints, lcp->bound->numvariablecharges);
  Matrix_transpose(&Actranspose, lcp->numconstraints, lcp->bound->numvariablecharges);

  Matrix_copy(Lhatlocal, lcp->Lhat, lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
  if (lcp->penaltyMatrix != NULL) { 
	 Matrix_add(Lhatlocal, Lhatlocal, lcp->penaltyMatrix, lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
  }
  
  Matrix_copypiece(LhatAndAc, 0, 0, Lhatlocal, 0, 0,
						 lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
  Matrix_copypiece(LhatAndAc, 0, lcp->bound->numvariablecharges, Actranspose, 0, 0,
						 lcp->bound->numvariablecharges, lcp->numconstraints);
  Matrix_copypiece(LhatAndAc, lcp->bound->numvariablecharges, 0, lcp->A_c, 0, 0,
						 lcp->numconstraints, lcp->bound->numvariablecharges);
  Matrix_pseudoinverse(lcp->LhatAndAcInv, LhatAndAc, lcp->bound->numvariablecharges + lcp->numconstraints,
							  lcp->bound->numvariablecharges + lcp->numconstraints);
  
  Matrix_free(Actranspose);
  Matrix_free(LhatAndAc);
  Matrix_free(Lhatlocal);
}

void LinConstrainedProblem_setupPreconditionerOrig(LinConstrainedProblem lcp) {
   Matrix Lhat, LhatAndAc,Actranspose;

   if (lcp->LhatAndAcInv != NULL)
      Matrix_free(lcp->LhatAndAcInv);
   
   lcp->LhatAndAcInv = Matrix_allocate(lcp->bound->numvariablecharges + lcp->numconstraints,
                                       lcp->bound->numvariablecharges + lcp->numconstraints);
   Lhat = Matrix_allocate(lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
   LhatAndAc         = Matrix_allocate(lcp->bound->numvariablecharges + lcp->numconstraints,
                                       lcp->bound->numvariablecharges + lcp->numconstraints);
   Actranspose = Matrix_allocate(lcp->numconstraints, lcp->bound->numvariablecharges);
   
   Matrix_copy(Actranspose, lcp->A_c, lcp->numconstraints, lcp->bound->numvariablecharges);
   Matrix_transpose(&Actranspose, lcp->numconstraints, lcp->bound->numvariablecharges);
   
   Optimizer_computePreconditionerHessian(lcp->bound, lcp->unbound, Lhat);
   Matrix_copypiece(LhatAndAc, 0, 0, Lhat, 0, 0,
                    lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
   Matrix_copypiece(LhatAndAc, 0, lcp->bound->numvariablecharges, Actranspose, 0, 0,
                    lcp->bound->numvariablecharges, lcp->numconstraints);
   Matrix_copypiece(LhatAndAc, lcp->bound->numvariablecharges, 0, lcp->A_c, 0, 0,
                    lcp->numconstraints, lcp->bound->numvariablecharges);
   Matrix_pseudoinverse(lcp->LhatAndAcInv, LhatAndAc, lcp->bound->numvariablecharges + lcp->numconstraints,
                        lcp->bound->numvariablecharges + lcp->numconstraints);
   
   Matrix_free(Actranspose);
   Matrix_free(LhatAndAc);
   Matrix_free(Lhat);
}

void LinConstrainedProblem_solve(LinConstrainedProblem lcp, Vector optimalCharges,
                                 Vector lagrangeMultipliers, Vector **approxCharges,
											int *numAdjustments, int* indexLooser, int* indexTighter) {
   unsigned int matrix_size = lcp->numconstraints + lcp->bound->numvariablecharges
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   unsigned int i;
   FILE *junk;
   Vector RHS = Vector_allocate(matrix_size);
   Vector soln = Vector_allocate(matrix_size);

   if (lcp->A_c == NULL || lcp->b == NULL) {
      printf("LinConstrainedProblem_solve: A_c matrix or b vector is null!  Exiting...\n");
      exit(-1);
   }

	Matrix Lhat;
	Lhat = Matrix_allocate(lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
	Optimizer_computePreconditionerHessian(lcp->bound, lcp->unbound, Lhat);
	if (lcp->LhatAndAcInv == NULL) // free it antime you want/need it reset, ie set penalty or change constraints
	  LinConstrainedProblem_setupPreconditioner(lcp);

   printf("setting up rhs\n");
   LinConstrainedProblem_setupRHS(lcp, RHS);

/*    junk = fopen("rhs","w"); */
/*    for (i = 0; i < matrix_size; i++) */
/*       fprintf(junk, "%f\n", RHS[i]); */
/*    fclose(junk); */
   printf("calling GMRES...\n");
   LinConstrainedProblem_GMRES(lcp, soln, RHS);
   
   Vector_copypiece(optimalCharges, 0, soln, 0, lcp->bound->numvariablecharges);
   Vector_copypiece(lagrangeMultipliers, 0, soln, lcp->bound->numvariablecharges, lcp->numconstraints);
   
   // clean up, go home
   Vector_free(RHS);
   Vector_free(soln);
	Matrix_free(Lhat);
}

void LinConstrainedProblem_setupRHS(LinConstrainedProblem lcp, Vector RHS) {
   unsigned int matrix_size = lcp->numconstraints + lcp->bound->numvariablecharges
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   Vector_zero(RHS, matrix_size);
   if (lcp->linearTerm == NULL) {
      lcp->linearTerm = Vector_allocate(lcp->bound->numvariablecharges);
      Optimizer_computeLinearTerm(lcp->bound, lcp->unbound, lcp->linearTerm);
   }
   // must put -c on RHS
   Vector_scale(lcp->linearTerm, -1.0, lcp->bound->numvariablecharges);
   Vector_copypiece(RHS, 0, lcp->linearTerm, 0, lcp->bound->numvariablecharges);
   Vector_scale(lcp->linearTerm, -1.0, lcp->bound->numvariablecharges);
   Vector_copypiece(RHS, lcp->bound->numvariablecharges, lcp->b, 0, lcp->numconstraints);
}

void LinConstrainedProblem_GMRES(LinConstrainedProblem lcp, Vector sol, Vector rhs) {
   unsigned int size = lcp->bound->numvariablecharges + lcp->numconstraints
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   real normr;
   unsigned int i;
   int j, k;
   real residual;

   r = Vector_allocate(size);

   LinConstrainedProblem_preconditionerMultiply(lcp, r, rhs);

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
      LinConstrainedProblem_operatorMultiply(lcp, P, bv[i]);
      LinConstrainedProblem_preconditionerMultiply(lcp, P, P);
      
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

void LinConstrainedProblem_operatorMultiply(LinConstrainedProblem lcp, Vector Ax, Vector x) {
   unsigned int matrix_size = lcp->numconstraints + lcp->bound->numvariablecharges
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   Vector q_temp, lagrange_temp, phi_b_temp, phi_u_temp;
   Vector Aq_temp, Alagrange_temp, Aphi_b_temp, Aphi_u_temp;
   Vector Aq_work, Alagrange_work, Aphi_b_work, Aphi_u_work;
/*    printf("entering\n"); */
   q_temp  = Vector_allocate(lcp->bound->numvariablecharges);
   Aq_temp = Vector_allocate(lcp->bound->numvariablecharges);
   Aq_work = Vector_allocate(lcp->bound->numvariablecharges);
   lagrange_temp   = Vector_allocate(lcp->numconstraints);
   Alagrange_temp  = Vector_allocate(lcp->numconstraints);
   Alagrange_work  = Vector_allocate(lcp->numconstraints);
   phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   Aphi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   Aphi_b_work = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
   Aphi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
   Aphi_u_work = Vector_allocate(lcp->unbound->numtotalsurfacevariables);

   Vector_zero(Aq_temp, lcp->bound->numvariablecharges);
   Vector_zero(Alagrange_temp, lcp->numconstraints);
   Vector_zero(Aphi_b_temp, lcp->bound->numtotalsurfacevariables);
   Vector_zero(Aphi_u_temp, lcp->unbound->numtotalsurfacevariables);

/*    printf("splitting.\n"); */
   LinConstrainedProblem_split(lcp, q_temp, lagrange_temp, phi_b_temp, phi_u_temp, x);

/*    printf("row 1\n"); */
   // do stuff here
   //   ROW 1
	if (lcp->penaltyMatrix != NULL) {
	  Matrix_multiplyvector(Aq_temp, lcp->penaltyMatrix, q_temp, lcp->bound->numvariablecharges, lcp->bound->numvariablecharges);
	}
   Matrix_multiplyvector_transpose(Aq_work, lcp->A_c, lagrange_temp, lcp->numconstraints, lcp->bound->numvariablecharges);
   Vector_addvector(Aq_temp, Aq_work, lcp->bound->numvariablecharges);
   PBEproblem_applyA3(lcp->bound, Aq_work, phi_b_temp);
   Vector_addvector(Aq_temp, Aq_work, lcp->bound->numvariablecharges);
   PBEproblem_applyA3(lcp->unbound, Aq_work, phi_u_temp);
   Vector_addscaledvector(Aq_temp, -1.0, Aq_work, lcp->bound->numvariablecharges); // subtracting Lb-Lu!
/*    printf("row 2\n"); */
   //   ROW 2
   Matrix_multiplyvector(Alagrange_temp, lcp->A_c, q_temp, lcp->numconstraints, lcp->bound->numvariablecharges);

/*    printf("row 3\n"); */
   //   ROW 3
   PBEproblem_applyA1(lcp->bound, Aphi_b_work, q_temp);
   Vector_addscaledvector(Aphi_b_temp, -1.0, Aphi_b_work, lcp->bound->numtotalsurfacevariables); // A1, A2 have opposite signs
   PBEproblem_applyA2(lcp->bound, Aphi_b_work, phi_b_temp);
   Vector_addvector(Aphi_b_temp, Aphi_b_work, lcp->bound->numtotalsurfacevariables);

/*    printf("row 4\n"); */
   //   ROW 4
   PBEproblem_applyA1(lcp->unbound, Aphi_u_work, q_temp);
   Vector_addscaledvector(Aphi_u_temp, -1.0, Aphi_u_work, lcp->unbound->numtotalsurfacevariables); // A1, A2 have opposite signs
   PBEproblem_applyA2(lcp->unbound, Aphi_u_work, phi_u_temp);
   Vector_addvector(Aphi_u_temp, Aphi_u_work, lcp->unbound->numtotalsurfacevariables);
   // end do stuff here

/*    printf("joining\n"); */
   LinConstrainedProblem_join(lcp, Ax, Aq_temp, Alagrange_temp, Aphi_b_temp, Aphi_u_temp);

   // clean up, go home
   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Aq_temp);
   Vector_free(Alagrange_temp);
   Vector_free(Aphi_b_temp);
   Vector_free(Aphi_u_temp);
   Vector_free(Aq_work);
   Vector_free(Alagrange_work);
   Vector_free(Aphi_b_work);
   Vector_free(Aphi_u_work);
/*    printf("returning\n"); */
}


void LinConstrainedProblem_preconditionerMultiply(LinConstrainedProblem lcp, Vector Px, Vector x) {
   Vector chargeAndLagrange, PchargeAndLagrange;
   Vector q_temp, lagrange_temp, phi_b_temp, phi_u_temp;
   Vector Pq_temp, Plagrange_temp, Pphi_b_temp, Pphi_u_temp;
   unsigned int sizeLhatandAc = lcp->bound->numvariablecharges + lcp->numconstraints;
   chargeAndLagrange = Vector_allocate(sizeLhatandAc);
   PchargeAndLagrange = Vector_allocate(sizeLhatandAc);
   q_temp  = Vector_allocate(lcp->bound->numvariablecharges);
   Pq_temp = Vector_allocate(lcp->bound->numvariablecharges);
   lagrange_temp   = Vector_allocate(lcp->numconstraints);
   Plagrange_temp  = Vector_allocate(lcp->numconstraints);
   phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   Pphi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
   Pphi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
	Vector P1chargeAndLagrange, P1phi_b_temp, P1phi_u_temp;
   P1chargeAndLagrange = Vector_allocate(sizeLhatandAc);
   P1phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   P1phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
	Vector P2P1chargeAndLagrange, P2P1phi_b_temp, P2P1phi_u_temp;
   P2P1chargeAndLagrange = Vector_allocate(sizeLhatandAc);
   P2P1phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   P2P1phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
	Vector P3P2P1chargeAndLagrange, P3P2P1phi_b_temp, P3P2P1phi_u_temp;
   P3P2P1chargeAndLagrange = Vector_allocate(sizeLhatandAc);
   P3P2P1phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   P3P2P1phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
	Vector A1bq_temp, A1uq_temp;
	A1bq_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
	A1uq_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
	Vector A3P1phi_b, A3P1phi_u;
	A3P1phi_b = Vector_allocate(lcp->bound->numvariablecharges);
	A3P1phi_u = Vector_allocate(lcp->bound->numvariablecharges);

   LinConstrainedProblem_split(lcp, q_temp, lagrange_temp, phi_b_temp, phi_u_temp, x);
	Vector_copypiece(chargeAndLagrange, 0, q_temp, 0, lcp->bound->numvariablecharges);
	Vector_copypiece(chargeAndLagrange, lcp->bound->numvariablecharges, lagrange_temp, 0, lcp->numconstraints);
	
	// P1 : [I 0; 0 D2^{-1}]
	Vector_copy(P1chargeAndLagrange, chargeAndLagrange, sizeLhatandAc);
	PBEproblem_applyPreconditioner(lcp->bound, P1phi_b_temp, phi_b_temp);
	PBEproblem_applyPreconditioner(lcp->unbound, P1phi_u_temp, phi_u_temp);
	
	// P2 : [I -A3; 0 I] (requires some finesse about A3)
	PBEproblem_applyA3(lcp->bound, A3P1phi_b, P1phi_b_temp);
	PBEproblem_applyA3(lcp->unbound, A3P1phi_u, P1phi_u_temp);
	Vector_scale(A3P1phi_b, -1., lcp->bound->numvariablecharges);
	Vector_addvector(A3P1phi_b, A3P1phi_u, lcp->bound->numvariablecharges);
	Vector_zero(P2P1chargeAndLagrange, sizeLhatandAc);
	Vector_copypiece(P2P1chargeAndLagrange, 0, A3P1phi_b, 0, lcp->bound->numvariablecharges);
	Vector_addvector(P2P1chargeAndLagrange, P1chargeAndLagrange, sizeLhatandAc);
	Vector_copy(P2P1phi_b_temp, P1phi_b_temp, lcp->bound->numtotalsurfacevariables);
	Vector_copy(P2P1phi_u_temp, P1phi_u_temp, lcp->unbound->numtotalsurfacevariables);

	// P3 : [LcAhatInv 0; 0 I]
	Matrix_multiplyvector(P3P2P1chargeAndLagrange, lcp->LhatAndAcInv, P2P1chargeAndLagrange, sizeLhatandAc, sizeLhatandAc);
	Vector_copy(P3P2P1phi_b_temp, P2P1phi_b_temp, lcp->bound->numtotalsurfacevariables);
	Vector_copy(P3P2P1phi_u_temp, P2P1phi_u_temp, lcp->unbound->numtotalsurfacevariables);

	// P4 : [I 0; D2^{-1}*A1 I]
	Vector_copy(PchargeAndLagrange, P3P2P1chargeAndLagrange, sizeLhatandAc);
	Vector_copypiece(q_temp, 0, P3P2P1chargeAndLagrange, 0, lcp->bound->numvariablecharges);
	PBEproblem_applyA1(lcp->bound, A1bq_temp, q_temp);
	PBEproblem_applyA1(lcp->unbound, A1uq_temp, q_temp);
	PBEproblem_applyPreconditioner(lcp->bound, Pphi_b_temp, A1bq_temp);
	PBEproblem_applyPreconditioner(lcp->unbound, Pphi_u_temp, A1uq_temp);
	Vector_addvector(Pphi_b_temp, P3P2P1phi_b_temp, lcp->bound->numtotalsurfacevariables);
	Vector_addvector(Pphi_u_temp, P3P2P1phi_u_temp, lcp->unbound->numtotalsurfacevariables); 
	
	Vector_copypiece(Pq_temp, 0, PchargeAndLagrange, 0, lcp->bound->numvariablecharges);
	Vector_copypiece(Plagrange_temp, 0, PchargeAndLagrange, lcp->bound->numvariablecharges, lcp->numconstraints);
   LinConstrainedProblem_join(lcp, Px, Pq_temp, Plagrange_temp, Pphi_b_temp, Pphi_u_temp);

	Vector_free(P1chargeAndLagrange); Vector_free(P1phi_b_temp); Vector_free(P1phi_u_temp);
	Vector_free(P2P1chargeAndLagrange); Vector_free(P2P1phi_b_temp); Vector_free(P2P1phi_u_temp);
	Vector_free(P3P2P1chargeAndLagrange); Vector_free(P3P2P1phi_b_temp); Vector_free(P3P2P1phi_u_temp);
	Vector_free(A1bq_temp); Vector_free(A1uq_temp);
	Vector_free(A3P1phi_b); Vector_free(A3P1phi_u);

   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Pq_temp);
   Vector_free(Plagrange_temp);
   Vector_free(Pphi_b_temp);
   Vector_free(Pphi_u_temp);
   Vector_free(chargeAndLagrange);
   Vector_free(PchargeAndLagrange);
}

void LinConstrainedProblem_preconditionerMultiplyOrig(LinConstrainedProblem lcp, Vector Px, Vector x) {
   Vector chargeAndLagrange, PchargeAndLagrange;
   Vector q_temp, lagrange_temp, phi_b_temp, phi_u_temp;
   Vector Pq_temp, Plagrange_temp, Pphi_b_temp, Pphi_u_temp;
   unsigned int sizeLhatandAc = lcp->bound->numvariablecharges + lcp->numconstraints;
/*    printf("entering precond\n"); */
   chargeAndLagrange = Vector_allocate(sizeLhatandAc);
   PchargeAndLagrange = Vector_allocate(sizeLhatandAc);
   q_temp  = Vector_allocate(lcp->bound->numvariablecharges);
   Pq_temp = Vector_allocate(lcp->bound->numvariablecharges);
   lagrange_temp   = Vector_allocate(lcp->numconstraints);
   Plagrange_temp  = Vector_allocate(lcp->numconstraints);
   phi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   Pphi_b_temp = Vector_allocate(lcp->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);
   Pphi_u_temp = Vector_allocate(lcp->unbound->numtotalsurfacevariables);

/*    printf("splitting.\n");          */
   LinConstrainedProblem_split(lcp, q_temp, lagrange_temp, phi_b_temp, phi_u_temp, x);
   Vector_copypiece(chargeAndLagrange, 0, q_temp, 0, lcp->bound->numvariablecharges);
   Vector_copypiece(chargeAndLagrange, lcp->bound->numvariablecharges, lagrange_temp, 0, lcp->numconstraints);
   
   //   Matrix_multiplyvector(PchargeAndLagrange, lcp->LhatAndAcInv, chargeAndLagrange, sizeLhatandAc, sizeLhatandAc);
/*    printf("preconditioning rows 1 and 2\n"); */
   Vector_copy(PchargeAndLagrange, chargeAndLagrange, sizeLhatandAc);
/*    printf("preconditioning row 3\n"); */
   Preconditioner_solve(Pphi_b_temp, lcp->bound->preconditioner, phi_b_temp);
/*    printf("preconditioning row 4\n"); */
   Preconditioner_solve(Pphi_u_temp, lcp->unbound->preconditioner, phi_u_temp);

   Vector_copypiece(Pq_temp, 0, PchargeAndLagrange, 0, lcp->bound->numvariablecharges);
   Vector_copypiece(Plagrange_temp, 0, PchargeAndLagrange, lcp->bound->numvariablecharges, lcp->numconstraints);

   LinConstrainedProblem_join(lcp, Px, Pq_temp, Plagrange_temp, Pphi_b_temp, Pphi_u_temp);

   Vector_free(q_temp);
   Vector_free(lagrange_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Pq_temp);
   Vector_free(Plagrange_temp);
   Vector_free(Pphi_b_temp);
   Vector_free(Pphi_u_temp);
   Vector_free(chargeAndLagrange);
   Vector_free(PchargeAndLagrange);
/*    printf("exiting precond\n"); */
}

void LinConstrainedProblem_operatorSave(LinConstrainedProblem lcp, char *filename) {
   unsigned int i,j;
   unsigned int matrix_size = lcp->bound->numvariablecharges + lcp->numconstraints
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   FILE *OUT = NULL;
   Vector x, Ax;
   x = Vector_allocate(matrix_size);
   Ax = Vector_allocate(matrix_size);

   OUT = fopen(filename, "w");
   if (OUT == NULL) {
      printf("LinConstrainedProblem_operatorSave: Could not open file %s.  Exiting.\n", filename);
      exit(-1);
   }

   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;

      LinConstrainedProblem_operatorMultiply(lcp, Ax, x);

      for (j = 0; j < matrix_size; j++)
         fprintf(OUT, "%f  ", Ax[j]);
      fprintf(OUT, "\n");
   }

   fclose(OUT);
   Vector_free(x);
   Vector_free(Ax);
}

void LinConstrainedProblem_preconditionerSave(LinConstrainedProblem lcp, char *filename) {
   unsigned int i,j;
   unsigned int matrix_size = lcp->bound->numvariablecharges + lcp->numconstraints
      + (lcp->bound->numtotalsurfacevariables + lcp->unbound->numtotalsurfacevariables);
   FILE *OUT = NULL;
   Vector x, Px;
   x = Vector_allocate(matrix_size);
   Px = Vector_allocate(matrix_size);

   OUT = fopen(filename, "w");
   if (OUT == NULL) {
      printf("LinConstrainedProblem_preconditionerSave: Could not open file %s.  Exiting.\n", filename);
      exit(-1);
   }

   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;

      LinConstrainedProblem_preconditionerMultiply(lcp, Px, x);

      for (j = 0; j < matrix_size; j++)
         fprintf(OUT, "%f  ", Px[j]);
      fprintf(OUT, "\n");
   }

   fclose(OUT);
   Vector_free(x);
   Vector_free(Px);
}

void LinConstrainedProblem_join(LinConstrainedProblem lcp, Vector dest, Vector q_temp, Vector lagrange_temp, Vector phi_b_temp, Vector phi_u_temp) {
   Vector_copypiece(dest, 0, q_temp, 0, lcp->bound->numvariablecharges);
   Vector_copypiece(dest, lcp->bound->numvariablecharges, lagrange_temp, 0, lcp->numconstraints);
   Vector_copypiece(dest, lcp->bound->numvariablecharges + lcp->numconstraints,
                    phi_b_temp, 0, lcp->bound->numtotalsurfacevariables);
   Vector_copypiece(dest, lcp->bound->numvariablecharges + lcp->numconstraints + lcp->bound->numtotalsurfacevariables,
                    phi_u_temp, 0, lcp->unbound->numtotalsurfacevariables);
}

void LinConstrainedProblem_split(LinConstrainedProblem lcp, Vector q_temp, Vector lagrange_temp, Vector phi_b_temp, Vector phi_u_temp, Vector src) {
   Vector_copypiece(q_temp, 0, src, 0, lcp->bound->numvariablecharges);
   Vector_copypiece(lagrange_temp, 0, src, lcp->bound->numvariablecharges, lcp->numconstraints);
   Vector_copypiece(phi_b_temp, 0, src, lcp->bound->numvariablecharges + lcp->numconstraints,
                    lcp->bound->numtotalsurfacevariables);
   Vector_copypiece(phi_u_temp, 0, src, lcp->bound->numvariablecharges + lcp->numconstraints + lcp->bound->numtotalsurfacevariables,
                    lcp->unbound->numtotalsurfacevariables);
}


void LinConstrainedProblem_loadConstraints(LinConstrainedProblem lcp,
                                           char *filename) {
}
