#include "Unconstrained.h"

UnconstrainedProblem UnconstrainedProblem_allocate() {
   UnconstrainedProblem up;
   up = (UnconstrainedProblem)calloc(1, sizeof(_UnconstrainedProblem));
   up->LhatInv = NULL;
	up->Lhat = NULL;
   up->linearTerm = NULL;
	up->penaltyMatrix = NULL;
	up->rightSingularVectors = NULL;
	up->leftSingularVectors = NULL;
	up->singularValues = NULL;
   return up;
}

void UnconstrainedProblem_free(UnconstrainedProblem up) {
  Matrix_free(up->Lhat);
   Matrix_free(up->LhatInv);
	Matrix_free(up->leftSingularVectors);
	Matrix_free(up->rightSingularVectors);
	Matrix_free(up->penaltyMatrix);
   Vector_free(up->linearTerm);
	Vector_free(up->singularValues);
   free(up);
}

void UnconstrainedProblem_setPBEproblems(UnconstrainedProblem up, PBEproblem bound, PBEproblem unbound) {
   up->bound = bound;
   up->unbound = unbound;
}

void UnconstrainedProblem_setPenalty(UnconstrainedProblem up, penaltyToleranceType penaltyType,
												 real maxEig, real penalty, Matrix penaltyMatrix,
												 Matrix leftSingularVectors, Matrix rightSingularVectors,
												 Vector singularValues,
												 real tolerance, real looserTol, real tighterTol) {
  if (up->penaltyMatrix != NULL) {
	 Matrix_free(up->penaltyMatrix);
  }
  
  if (up->leftSingularVectors != NULL) {
	 Matrix_free(up->leftSingularVectors);
  }

  if (up->rightSingularVectors != NULL) {
	 Matrix_free(up->rightSingularVectors);
  }

  if (up->singularValues != NULL) {
	 Vector_free(up->singularValues);
  }
  
  up->penalty = penalty;
  up->penaltyMatrix = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  Matrix_copy(up->penaltyMatrix, penaltyMatrix, up->bound->numvariablecharges, up->bound->numvariablecharges);
  up->rightSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  Matrix_copy(up->rightSingularVectors, rightSingularVectors, up->bound->numvariablecharges, up->bound->numvariablecharges);
  up->leftSingularVectors = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  Matrix_copy(up->leftSingularVectors, leftSingularVectors, up->bound->numvariablecharges, up->bound->numvariablecharges);
  up->singularValues = Vector_allocate(up->bound->numvariablecharges);
  Vector_copy(up->singularValues, singularValues, up->bound->numvariablecharges);
  up->tolerance = tolerance;
  up->looserTol = looserTol;
  up->tighterTol = tighterTol;
  up->penaltyType = penaltyType;
  up->maxEigenvalue = maxEig;
}
  
void UnconstrainedProblem_setupPreconditioner(UnconstrainedProblem up) {
  Matrix Lhattmp = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  if (up->LhatInv != NULL) {  // will re-set at every solve.  maybe a little slow but at least guaranteed to do the right thing.
	  Matrix_free(up->LhatInv);
	  up->LhatInv = NULL;
  }
  
   up->LhatInv = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);

	if (up->Lhat == NULL) {
	  up->Lhat        = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
	  Optimizer_computePreconditionerHessian(up->bound, up->unbound, up->Lhat);
	}
	Matrix_copy(Lhattmp, up->Lhat, up->bound->numvariablecharges, up->bound->numvariablecharges);
	
	if (up->penaltyMatrix != NULL) {
	  Matrix_add(Lhattmp, up->penaltyMatrix, up->Lhat, up->bound->numvariablecharges, up->bound->numvariablecharges);
	}

   Matrix_pseudoinverse(up->LhatInv, Lhattmp, up->bound->numvariablecharges, up->bound->numvariablecharges);
	Matrix_free(Lhattmp);
}

void UnconstrainedProblem_solve(UnconstrainedProblem up, Vector optimalCharges,
										  Vector **approxCharges, int *numAdjustments,
										  int *lowerAdjustment, int* higherAdjustment) {
   unsigned int matrix_size = up->bound->numtotalsurfacevariables
      + up->unbound->numtotalsurfacevariables
      + up->bound->numvariablecharges;
   Vector RHS = Vector_allocate(matrix_size);
   Vector soln = Vector_allocate(matrix_size);

	UnconstrainedProblem_setupPreconditioner(up); // GETS CALLED WITH EVERY CALL TO _SOLVE!
   UnconstrainedProblem_setupRHS(up, RHS);
   
	FILE *rhs = fopen("rhs.m", "w");
	unsigned int i;
	for(i = 0; i < matrix_size; i++) {
	  fprintf(rhs, "%f\n", RHS[i]);
	}
	fclose(rhs);
   UnconstrainedProblem_GMRES(up, soln, RHS);

   Vector_copy(optimalCharges, soln, up->bound->numvariablecharges);
	unsigned int looserTolIndex, origTolIndex, tighterTolIndex;
	looserTolIndex = origTolIndex = tighterTolIndex = 0;
	if ((up->looserTol > up->tolerance) || (up->tighterTol < up->tolerance)) {
	  for (i = 0 ; i < up->bound->numvariablecharges; i++) {
		 if (up->singularValues[i] > up->looserTol * up->singularValues[0])
			looserTolIndex = i;
		 if (up->singularValues[i] > up->tolerance * up->singularValues[0])
			origTolIndex = i;
		 if (up->singularValues[i] > up->tighterTol * up->singularValues[0])
			tighterTolIndex = i;
	  }
	  looserTolIndex += 1; origTolIndex += 1; tighterTolIndex += 1;
	  printf("looserTolIndex = %d; origTolIndex = %d; tighterTolIndex = %d\n",
				looserTolIndex, origTolIndex, tighterTolIndex);
	  *lowerAdjustment = looserTolIndex;
	  *higherAdjustment = tighterTolIndex;
	  *numAdjustments = tighterTolIndex - looserTolIndex + 1;
	  *approxCharges = (Vector *)calloc(*numAdjustments, sizeof(Vector));
	  unsigned int curApproxIndex = 0;
	  for (i = looserTolIndex; i <= tighterTolIndex; i++) {
		 (*approxCharges)[curApproxIndex] = Vector_allocate(up->bound->numvariablecharges);
		 printf("doing problem unc_approx_%d\n", up->bound->numvariablecharges - i);
		 if (i-1 == origTolIndex) {
			printf("curApproxIndex = %d, i = %d\n", curApproxIndex, i);
		 }
		 UnconstrainedProblem_getApproxSWMSoln(up, (*approxCharges)[curApproxIndex], optimalCharges,
														 	origTolIndex, i);
		 curApproxIndex++;
	  }
	  Matrix junk = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
	  unsigned int j;
	  for (i = 0; i < up->bound->numvariablecharges; i++) {
		 Vector x = Vector_allocate(up->bound->numvariablecharges);
		 Vector Ax = Vector_allocate(up->bound->numvariablecharges);
		 Vector_zero(x, up->bound->numvariablecharges);
		 x[i] = 1.0;
		 UnconstrainedProblem_getApproxSWMSoln(up, Ax, x, origTolIndex, i+1);
		 for (j = 0; j < up->bound->numvariablecharges; j++)
			junk[j][i] = Ax[j];
		 Vector_free(x);
		 Vector_free(Ax);
	  }
	  Matrix_writefile("checkMatc", junk, up->bound->numvariablecharges, up->bound->numvariablecharges);
	  Matrix_free(junk);
	}
	
	
   // clean up, go home
   Vector_free(soln);
   Vector_free(RHS);
}

void UnconstrainedProblem_getApproxSWMSoln(UnconstrainedProblem up, Vector newsoln, Vector oldsoln,
														 unsigned int orig, unsigned int newInd) {
  //  newInd = orig+1;
  if (orig == newInd) {
	 Vector_copy(newsoln, oldsoln, up->bound->numvariablecharges);
	 return;
  }
  unsigned int i;
  unsigned int rankupdate = abs(orig - newInd);
  real penalty = up->penalty;
  Matrix mod = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  Matrix small = Matrix_allocate(rankupdate, rankupdate);
  Matrix Vk = Matrix_allocate(up->bound->numvariablecharges, rankupdate);
  Matrix Uk = Matrix_allocate(up->bound->numvariablecharges, rankupdate);
  Matrix VkT = Matrix_allocate(up->bound->numvariablecharges, rankupdate);
  Matrix Ainv_Uk = Matrix_allocate(up->bound->numvariablecharges, rankupdate);
  Matrix VkT_Ainv_Uk = Matrix_allocate(rankupdate, rankupdate);
  Matrix smallInv_VkT = Matrix_allocate(rankupdate, up->bound->numvariablecharges);
  Matrix Ainv_Uk_smallInv_VkT = Matrix_allocate(up->bound->numvariablecharges, up->bound->numvariablecharges);
  if (newInd < orig) {
	 Matrix_copypiece(Vk, 0, 0, up->rightSingularVectors, 0, newInd, up->bound->numvariablecharges, rankupdate);
	 Matrix_copypiece(Uk, 0, 0, up->leftSingularVectors, 0, newInd, up->bound->numvariablecharges, rankupdate);
  } else {
	 Matrix_copypiece(Vk, 0, 0, up->rightSingularVectors, 0, orig, up->bound->numvariablecharges, rankupdate);
	 Matrix_copypiece(Uk, 0, 0, up->leftSingularVectors, 0, orig, up->bound->numvariablecharges, rankupdate);
	 penalty = - penalty;
  }
  Matrix_copy(VkT, Vk, up->bound->numvariablecharges, rankupdate);
  Matrix_transpose(&VkT, up->bound->numvariablecharges, rankupdate);
  Matrix_free(Vk);

  Matrix_scale(Uk, penalty, rankupdate, up->bound->numvariablecharges);

  Matrix_multiplymatrix(Ainv_Uk, up->LhatInv, Uk, up->bound->numvariablecharges, up->bound->numvariablecharges, rankupdate);
  Matrix_writefile("LhatInv", up->LhatInv, up->bound->numvariablecharges, up->bound->numvariablecharges);
  Matrix_writefile("AinvUk", Ainv_Uk, up->bound->numvariablecharges, rankupdate);
  Matrix_multiplymatrix(VkT_Ainv_Uk, VkT, Ainv_Uk, rankupdate, up->bound->numvariablecharges, rankupdate);
  for (i = 0; i < rankupdate; i++)
	 VkT_Ainv_Uk[i][i] += 1.0;  // identity plus
  Matrix_pseudoinverse(small, VkT_Ainv_Uk, rankupdate, rankupdate);

  Matrix_multiplymatrix(smallInv_VkT, small, VkT, rankupdate, rankupdate, up->bound->numvariablecharges);
  Matrix_multiplymatrix(Ainv_Uk_smallInv_VkT, Ainv_Uk, smallInv_VkT, up->bound->numvariablecharges, rankupdate, up->bound->numvariablecharges);
  Matrix_scale(Ainv_Uk_smallInv_VkT, -1.0, up->bound->numvariablecharges, up->bound->numvariablecharges);
  for (i = 0; i < up->bound->numvariablecharges; i++)
	 Ainv_Uk_smallInv_VkT[i][i] += 1.0;
  Vector_zero(newsoln, up->bound->numvariablecharges);
/*   Matrix_writefile("Ainv", Ainv_Uk_smallInv_VkT, up->bound->numvariablecharges, up->bound->numvariablecharges); */
  Matrix_multiplyvector(newsoln, Ainv_Uk_smallInv_VkT, oldsoln, up->bound->numvariablecharges, up->bound->numvariablecharges);
/*   Vector_writefile("oldsoln", oldsoln, up->bound->numvariablecharges); */
/*   Vector_writefile("newsoln", newsoln, up->bound->numvariablecharges); */
  // clean up, go home
  Matrix_free(mod);
  Matrix_free(small);
  Matrix_free(Uk);
  Matrix_free(VkT);
  Matrix_free(Ainv_Uk);
  Matrix_free(VkT_Ainv_Uk);
  Matrix_free(smallInv_VkT);
  Matrix_free(Ainv_Uk_smallInv_VkT);
}

void UnconstrainedProblem_GMRES(UnconstrainedProblem up, Vector sol, Vector rhs) {
   unsigned int size = up->bound->numvariablecharges +
      + (up->bound->numtotalsurfacevariables 
         + up->unbound->numtotalsurfacevariables);
   Vector r, x, c, s, g, y, P, bv[MAXITERATIONS+1];
   Matrix H;
   Matrix Y, BV, Intermed;
   real normr;
   unsigned int i;
   int j, k;
   real residual;
	printf("in UnconstrainedProblem_GMRES!\n");
   r = Vector_allocate(size);

   UnconstrainedProblem_preconditionerMultiply(up, r, rhs);

   normr = Vector_norm(r, size);
	printf("normr = %f\n", normr);
	
   x = Vector_allocate(size);

   c = Vector_allocate(MAXITERATIONS+1);
   s = Vector_allocate(MAXITERATIONS+1);
   g = Vector_allocate(MAXITERATIONS+1);
   y = Vector_allocate(MAXITERATIONS+1);
   H = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   Y = Matrix_allocate(MAXITERATIONS+1, MAXITERATIONS+1);
   Intermed = Matrix_allocate(MAXITERATIONS+1, size);
   
   P = Vector_allocate(size);

   g[0] = Vector_norm(r, size);
   bv[0] = Vector_allocate(size);
   Vector_copy(bv[0], r, size);
   Vector_scale(bv[0], 1.0 / g[0], size);
   Vector_free(r);
   for (i = 0; i < MAXITERATIONS; i++) {
      UnconstrainedProblem_operatorMultiply(up, P, bv[i]);
      UnconstrainedProblem_preconditionerMultiply(up, P, P);
      
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

      for (k = 0; k < i; k++) // "looking inside GMRES" so we can recover all iterates later
         y[k] = g[k];

      
      g[i+1] = 0.0;
      givensrotate(c[i], s[i], &g[i], &g[i+1]);
      
      residual = fabs(g[i+1]) / normr;

#ifdef PRINT_GMRES_RESIDUALS
      printf("Iteration: %u Residual: %2.8f\n", i+1, residual);
#endif

      // begin save intermediate
		if (saveGMRES) {
      for (k = 0; k <= i; k++) // "looking inside GMRES" so we can recover all iterates later
         Y[k][i] = g[k];

      for (k = i; k >= 0; k--) {
         Y[k][i] /= H[k][k];
         for (j = k-1; j >= 0; j--)
            Y[j][i] -= H[k][j] * Y[k][i];
      }
      
      for (j = 0; j <= i; j++)
         Vector_addscaledvector(Intermed[i], Y[j][i], bv[j], size);
      }
      // end save intermediate

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

   // begin save intermediate stuff
	if (saveGMRES) {
   Matrix_writefile("Y", Y, i+1, i+1);
   BV = Matrix_allocate(i+1,size);
   for (j = 0; j <= i; j++)
      Vector_copy(BV[j], bv[j], size);
	Matrix_writefile("BV", BV, i+1, size);
   Matrix_writefile("H", H, i+1, i+2);
   Matrix_writefile("X", Intermed, i+1, up->bound->numvariablecharges);
	// the above is because we don't need all the rows of X, just the ones corresponding to the charge variables
	// otherwise X is an ENORMOUS file!
	}
   // end save intermediate stuff
   
   Vector_copy(sol, x, size);

   Vector_free(x);
   Vector_free(P);
   Vector_free(c);
   Vector_free(s);
   Vector_free(g);
   Vector_free(y);
   Matrix_free(H);
   
   for (j = 0; j <= i+1; j++)
      Vector_free(bv[j]);
}

void UnconstrainedProblem_setupRHS(UnconstrainedProblem up, Vector RHS) {
   if (up->linearTerm == NULL) {
       up->linearTerm = Vector_allocate(up->bound->numvariablecharges);
       Optimizer_computeLinearTerm(up->bound, up->unbound, up->linearTerm);
   }

   Vector_zero(RHS, (up->bound->numtotalsurfacevariables + up->unbound->numtotalsurfacevariables) + up->bound->numvariablecharges);
   Vector_scale(up->linearTerm, -1.0, up->bound->numvariablecharges);
   Vector_copypiece(RHS, 0, up->linearTerm, 0, up->bound->numvariablecharges);
   Vector_scale(up->linearTerm, -1.0, up->bound->numvariablecharges);
}

void UnconstrainedProblem_loadNewChargeDistribution(UnconstrainedProblem up, char *PDBfilename, char *CRGfilename) {
  PBEproblem_loadNewChargeDistribution(up->bound, PDBfilename, CRGfilename);
  PBEproblem_loadNewChargeDistribution(up->unbound, PDBfilename, CRGfilename);
  if (up->linearTerm != NULL) {
	 Vector_free(up->linearTerm);
	 up->linearTerm = Vector_allocate(up->bound->numvariablecharges);
  }
  if (up->LhatInv != NULL) {
	 Matrix_free(up->LhatInv);
  }
}


void UnconstrainedProblem_operatorMultiply(UnconstrainedProblem up, Vector Ax, Vector x) {
   Vector q_temp, phi_b_temp, phi_u_temp;
   Vector Aq_temp,  Aphi_b_temp,  Aphi_u_temp;
   Vector Aq_work,  Aphi_b_work,  Aphi_u_work;
   q_temp = Vector_allocate(up->bound->numvariablecharges);
   phi_b_temp = Vector_allocate(up->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(up->unbound->numtotalsurfacevariables);
   Aq_temp = Vector_allocate(up->bound->numvariablecharges);
   Aphi_b_temp = Vector_allocate(up->bound->numtotalsurfacevariables);
   Aphi_u_temp = Vector_allocate(up->unbound->numtotalsurfacevariables);
   Aq_work = Vector_allocate(up->bound->numvariablecharges);
   Aphi_b_work = Vector_allocate(up->bound->numtotalsurfacevariables);
   Aphi_u_work = Vector_allocate(up->unbound->numtotalsurfacevariables);

   Vector_zero(Aq_temp, up->bound->numvariablecharges);
   Vector_zero(Aphi_b_temp, up->bound->numtotalsurfacevariables);
   Vector_zero(Aphi_u_temp, up->unbound->numtotalsurfacevariables);
   
   UnconstrainedProblem_split(up, q_temp, phi_b_temp, phi_u_temp, x);

   // START MULT STUFF
   
   // row 1
	if (up->penaltyMatrix != NULL) {
	  Matrix_multiplyvector(Aq_work, up->penaltyMatrix, q_temp, up->bound->numvariablecharges, up->bound->numvariablecharges);
	  Vector_addvector(Aq_temp, Aq_work, up->bound->numvariablecharges);
	}
	//	printf("applying A3 operators\n");
   PBEproblem_applyA3(up->bound, Aq_work, phi_b_temp);
   Vector_addvector(Aq_temp, Aq_work, up->bound->numvariablecharges);
   PBEproblem_applyA3(up->unbound, Aq_work, phi_u_temp);
   Vector_addscaledvector(Aq_temp, -1.0, Aq_work, up->bound->numvariablecharges); // subtracting Lb - Lu!!

   // row 2: remember A1, A2 have opposite signs because we want A2 phi = A1 q
	//	printf("applying bound A1, A2\n");
   PBEproblem_applyA1(up->bound, Aphi_b_work, q_temp);
   Vector_addscaledvector(Aphi_b_temp, -1.0, Aphi_b_work, up->bound->numtotalsurfacevariables);
   PBEproblem_applyA2(up->bound, Aphi_b_work, phi_b_temp);
   Vector_addvector(Aphi_b_temp, Aphi_b_work, up->bound->numtotalsurfacevariables);

   // row 3: A1, A2 have opposite signs again of course
	//	printf("applying unbound A1, A2\n");
   PBEproblem_applyA1(up->unbound, Aphi_u_work, q_temp);
   Vector_addscaledvector(Aphi_u_temp, -1.0, Aphi_u_work, up->unbound->numtotalsurfacevariables);
   PBEproblem_applyA2(up->unbound, Aphi_u_work, phi_u_temp);
   Vector_addvector(Aphi_u_temp, Aphi_u_work, up->unbound->numtotalsurfacevariables);
   
   // END MULT STUFF
   
   UnconstrainedProblem_join(up, Ax, Aq_temp, Aphi_b_temp, Aphi_u_temp);

   Vector_free(q_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Aq_temp);
   Vector_free(Aphi_b_temp);
   Vector_free(Aphi_u_temp);
   Vector_free(Aq_work);
   Vector_free(Aphi_b_work);
   Vector_free(Aphi_u_work);
}

void UnconstrainedProblem_preconditionerMultiply(UnconstrainedProblem up, Vector Px, Vector x) {

  Vector q, phi_b, phi_u;
  Vector phi_b_temp, phi_u_temp, P_phi_b, P_phi_u;

  // identity preconditioner: should be commented out most of the time.
/*   Vector_copy(Px, x, up->bound->numvariablecharges + up->bound->numtotalsurfacevariables + */
/* 				  up->unbound->numtotalsurfacevariables); */
/*   return; */
  // end identity preconditioner
  
  q = Vector_allocate(up->bound->numvariablecharges);
  phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  phi_b_temp = Vector_allocate(up->bound->numtotalsurfacevariables);
  phi_u_temp = Vector_allocate(up->unbound->numtotalsurfacevariables);
  P_phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  P_phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  UnconstrainedProblem_preconditionerMultiplyLowerTriang(up, Px, x);  // does the first three products

  // return here for P3P2P1 (or less!) preconditioners (triangularize at most)
  //   return;
  
  // now multiply by [I 0; A2^{-1}*A1 I];
  UnconstrainedProblem_split(up, q, phi_b, phi_u, Px);

  PBEproblem_applyA1(up->bound, phi_b_temp, q);
  PBEproblem_applyA1(up->unbound, phi_u_temp, q);
  Preconditioner_solve(P_phi_b, up->bound->preconditioner, phi_b_temp);
  Preconditioner_solve(P_phi_u, up->unbound->preconditioner, phi_u_temp);
  Vector_addvector(P_phi_b, phi_b, up->bound->numtotalsurfacevariables);
  Vector_addvector(P_phi_u, phi_u, up->unbound->numtotalsurfacevariables);

  // now clean up, go home
  UnconstrainedProblem_join(up, Px, q, P_phi_b, P_phi_u);
  Vector_free(q);
  Vector_free(phi_b);
  Vector_free(phi_u);
  Vector_free(phi_b_temp);
  Vector_free(phi_u_temp);
  Vector_free(P_phi_b);
  Vector_free(P_phi_u);
}

void UnconstrainedProblem_preconditionerMultiplyLowerTriang(UnconstrainedProblem up, Vector Px, Vector x) {
  Vector q, phi_u, phi_b;
  Vector P1_q, P1_phi_u, P1_phi_b;
  Vector P2P1_q, P2P1_phi_u, P2P1_phi_b;
  Vector P_q, P_phi_u, P_phi_b;

  q = Vector_allocate(up->bound->numvariablecharges);
  P1_q = Vector_allocate(up->bound->numvariablecharges);
  P2P1_q = Vector_allocate(up->bound->numvariablecharges);
  P_q = Vector_allocate(up->bound->numvariablecharges);
  phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  P1_phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  P2P1_phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  P_phi_b = Vector_allocate(up->bound->numtotalsurfacevariables);
  phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  P1_phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  P2P1_phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  P_phi_u = Vector_allocate(up->unbound->numtotalsurfacevariables);
  
  UnconstrainedProblem_split(up, q, phi_b, phi_u, x);

  // see notes from 2/20-2/23 on "reduced space preconditioning"
  Vector_copy(P1_q, q, up->bound->numvariablecharges);
  Preconditioner_solve(P1_phi_b, up->bound->preconditioner, phi_b);
  Preconditioner_solve(P1_phi_u, up->unbound->preconditioner, phi_u);

  Vector_copy(P2P1_q, P1_q, up->bound->numvariablecharges);
  PBEproblem_applyA3(up->bound, q, P1_phi_b);
  PBEproblem_applyA3(up->unbound, P_q, P1_phi_u);
  Vector_scale(q,-1,up->bound->numvariablecharges);
  Vector_addvector(P2P1_q, q, up->bound->numvariablecharges);
  Vector_addvector(P2P1_q, P_q, up->bound->numvariablecharges);
  Vector_copy(P2P1_phi_b, P1_phi_b, up->bound->numtotalsurfacevariables);
  Vector_copy(P2P1_phi_u, P1_phi_u, up->unbound->numtotalsurfacevariables);

  Vector_zero(P_q, up->bound->numvariablecharges);
  Matrix_multiplyvector(P_q, up->LhatInv,  P2P1_q, up->bound->numvariablecharges, up->bound->numvariablecharges);
  Vector_copy(P_phi_u, P2P1_phi_u, up->unbound->numtotalsurfacevariables);
  Vector_copy(P_phi_b, P2P1_phi_b, up->bound->numtotalsurfacevariables);
    
  UnconstrainedProblem_join(up, Px, P_q, P_phi_b, P_phi_u);

  Vector_free(q); Vector_free(phi_u); Vector_free(phi_b);
  Vector_free(P1_q); Vector_free(P1_phi_u); Vector_free(P1_phi_b);
  Vector_free(P2P1_q); Vector_free(P2P1_phi_u); Vector_free(P2P1_phi_b);
  Vector_free(P_q); Vector_free(P_phi_u); Vector_free(P_phi_b);
}

void UnconstrainedProblem_preconditionerMultiplyOrig(UnconstrainedProblem up, Vector Px, Vector x) {
   Vector q_temp, Pq_temp, phi_b_temp, Pphi_b_temp, phi_u_temp, Pphi_u_temp;
   q_temp = Vector_allocate(up->bound->numvariablecharges);
   Pq_temp = Vector_allocate(up->bound->numvariablecharges);
   phi_b_temp = Vector_allocate(up->bound->numtotalsurfacevariables);
   Pphi_b_temp = Vector_allocate(up->bound->numtotalsurfacevariables);
   phi_u_temp = Vector_allocate(up->unbound->numtotalsurfacevariables);
   Pphi_u_temp = Vector_allocate(up->unbound->numtotalsurfacevariables);

   UnconstrainedProblem_split(up, q_temp, phi_b_temp, phi_u_temp, x);
   // do mult stuff here
   // row 1:
/*    Matrix_multiplyvector(Pq_temp, up->LhatInv, q_temp, up->bound->numvariablecharges, */
/*                          up->bound->numvariablecharges); */
   Vector_copy(Pq_temp, q_temp, up->bound->numvariablecharges);
   // row 2:
   Preconditioner_solve(Pphi_b_temp, up->bound->preconditioner, phi_b_temp);
   // row 3:
   Preconditioner_solve(Pphi_u_temp, up->unbound->preconditioner, phi_u_temp);
   // end mult stuff
   UnconstrainedProblem_join(up, Px, Pq_temp, Pphi_b_temp, Pphi_u_temp);

   Vector_free(q_temp);
   Vector_free(phi_b_temp);
   Vector_free(phi_u_temp);
   Vector_free(Pq_temp);
   Vector_free(Pphi_b_temp);
   Vector_free(Pphi_u_temp);
}

void UnconstrainedProblem_operatorSave(UnconstrainedProblem up, char *filename) {
   unsigned int i,j;
   unsigned int matrix_size = up->bound->numvariablecharges
      + (up->bound->numtotalsurfacevariables + up->unbound->numtotalsurfacevariables);
   FILE *OUT = NULL;
   Vector x, Ax;
   x = Vector_allocate(matrix_size);
   Ax = Vector_allocate(matrix_size);

   OUT = fopen(filename, "w");
   if (OUT == NULL) {
      printf("UnconstrainedProblem_operatorSave: Could not open file %s.  Exiting.\n", filename);
      exit(-1);
   }

   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;

      UnconstrainedProblem_operatorMultiply(up, Ax, x);

      for (j = 0; j < matrix_size; j++)
         fprintf(OUT, "%f  ", Ax[j]);
      fprintf(OUT, "\n");
   }

   fclose(OUT);
   Vector_free(x);
   Vector_free(Ax);
}

void UnconstrainedProblem_preconditionerSave(UnconstrainedProblem up, char *filename) {
   unsigned int i,j;
   unsigned int matrix_size = up->bound->numvariablecharges
      + (up->bound->numtotalsurfacevariables + up->unbound->numtotalsurfacevariables);
   FILE *OUT = NULL;
   Vector x, Px;
   x = Vector_allocate(matrix_size);
   Px = Vector_allocate(matrix_size);

   OUT = fopen(filename, "w");
   if (OUT == NULL) {
      printf("UnconstrainedProblem_operatorSave: Could not open file %s.  Exiting.\n", filename);
      exit(-1);
   }

   for (i = 0; i < matrix_size; i++) {
      Vector_zero(x, matrix_size);
      x[i] = 1.0;

      UnconstrainedProblem_preconditionerMultiply(up, Px, x);

      for (j = 0; j < matrix_size; j++)
         fprintf(OUT, "%f  ", Px[j]);
      fprintf(OUT, "\n");
   }

   fclose(OUT);
   Vector_free(x);
   Vector_free(Px);
}

void UnconstrainedProblem_join(UnconstrainedProblem up,
                        Vector dest, Vector q_temp, Vector phi_b_temp, Vector phi_u_temp) {
   Vector_copypiece(dest, 0, q_temp, 0, up->bound->numvariablecharges);
   Vector_copypiece(dest, up->bound->numvariablecharges, phi_b_temp, 0, up->bound->numtotalsurfacevariables);
   Vector_copypiece(dest, up->bound->numvariablecharges + up->bound->numtotalsurfacevariables,
                    phi_u_temp, 0, up->unbound->numtotalsurfacevariables);
}
 
void UnconstrainedProblem_split(UnconstrainedProblem up,
                         Vector q_temp, Vector phi_b_temp, Vector phi_u_temp, Vector src) {
   Vector_copypiece(q_temp, 0, src, 0, up->bound->numvariablecharges);
   Vector_copypiece(phi_b_temp, 0, src, up->bound->numvariablecharges, up->bound->numtotalsurfacevariables);
   Vector_copypiece(phi_u_temp, 0, src, up->bound->numvariablecharges + up->bound->numtotalsurfacevariables,
                    up->unbound->numtotalsurfacevariables);
}
