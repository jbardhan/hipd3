#include "Optimizer.h"
#include <sys/time.h>
#include <sys/resource.h>

extern char variablechain;

void dsyevr_(char* JOBZ, char* RANGE, char* UPLO, unsigned int* N,
				 double* A, unsigned int* LDA, double* VL, double* VU,
				 unsigned int* IL, unsigned int* IU, double* ABSTOL,
				 unsigned int* M, double* W, double* Z, unsigned int* LDZ,
				 unsigned int* ISUPPZ, double* WORK, unsigned int* LWORK,
				 unsigned int* IWORK, unsigned int* LIWORK, int* INFO);
void ssyevr_(char* JOBZ, char* RANGE, char* UPLO, unsigned int* N,
				 float* A, unsigned int* LDA, float* VL, float* VU,
				 unsigned int* IL, unsigned int* IU, float* ABSTOL,
				 unsigned int* M, float* W, float* Z, unsigned int* LDZ,
				 unsigned int* ISUPPZ, float* WORK, unsigned int* LWORK,
				 unsigned int* IWORK, unsigned int* LIWORK, int* INFO);

void Optimizer_computeHessian(PBEproblem bound, PBEproblem unbound, Matrix L,
										unsigned int startcol, unsigned int endcol) {
  unsigned int columnCount;
  unsigned int n_c = bound->numvariablecharges;
  unsigned int i;
  Vector chargeVec = Vector_allocate(n_c);
  Vector LchargeVec = Vector_allocate(n_c);
  unsigned int total_GMRES_iter = 0;

  struct rusage ruse;
  struct timeval tval;
  real starttime, endtime, startwalltime, endwalltime;

   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

  for (columnCount = startcol; columnCount < ((endcol>=n_c)?n_c:endcol+1);
		 columnCount++) { 
	 Vector_zero(chargeVec, n_c);
	 chargeVec[columnCount] = 1.0;
	 
	 Optimizer_multiplyByL(bound, unbound, LchargeVec, chargeVec);
	 total_GMRES_iter += num_GMRES_iter;

	 for (i = 0; i < n_c; i++) {
		L[i][columnCount] = LchargeVec[i];
	 }
	 Vector_zero(LchargeVec, n_c);
  }
  num_GMRES_iter = total_GMRES_iter;
  printf("total number of GMRES iterations required for Hessian: %d\n",
			num_GMRES_iter);

   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
   printf("time to compute Hessian: %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);

  Vector_free(chargeVec);
  Vector_free(LchargeVec);
}

void Optimizer_computePenaltyMatrix(PBEproblem bound, PBEproblem unbound, unsigned int problemsize, Matrix Lhat, 
												penaltyToleranceType penaltyType, real tolerance, real penaltyScale, Matrix penalty,
												Matrix leftSingularVectors, Matrix rightSingularVectors, Vector singularValues, real *maxEig) {
  // this function is basically Matrix_pseudoinverse_droptol, except
  // the SMALLEST singular vectors are retained--to be penalized in
  //  printf("in cpm tolerance = %f\n");
  char JOBZ = 'V', RANGE = 'A', UPLO = 'U';
  unsigned int N = problemsize;
  Matrix A = Matrix_allocate(N,N);
  unsigned int LDA = N, LDZ = N;
  real VL,VU;
  unsigned int IL, IU;
  real ABSTOL = -1.;
  unsigned int M;
  Matrix S = Matrix_allocate(N, N);
  Matrix V = Matrix_allocate(N,N);
  Matrix Vt = Matrix_allocate(N,N);
  Matrix temp = Matrix_allocate(N,N);
  Vector D = Vector_allocate(N);
  unsigned int *ISUPPZ = (unsigned int *)malloc(N * sizeof(unsigned int));
  real *WORK = (real *)malloc(26 * N * sizeof(real));
  unsigned int LWORK = 26 * N;
  unsigned int *IWORK = (unsigned int *)malloc(10 * N * sizeof(unsigned int));
  unsigned int LIWORK = 10 * N;
  int INFO;

  unsigned int i, j;
  for (i = 0; i < N; i++)
	 for (j = 0; j < N; j++)
		A[i][j] = .5 * (Lhat[i][j] + Lhat[j][i]);

#ifdef REAL_IS_DOUBLE
#ifdef OMP
#pragma omp critical (lapack)
#endif
/*   printf("JOBZ = %c\nRANGE = %c\nUPLO =  %c\nN = %d\nA = %d\nLDA = %d\nVL = %f\nVU = %f\n", JOBZ, RANGE, UPLO, N, A[0], LDA, VL, VU); */
  dsyevr_(&JOBZ, &RANGE, &UPLO, &N,
			 A[0], &LDA, &VL, &VU,
			 &IL, &IU, &ABSTOL,
			 &M, D, Vt[0], &LDZ,
			 ISUPPZ, WORK, &LWORK,
			 IWORK, &LIWORK, &INFO);  // V is returned as V'!!!
#else
#ifdef REAL_IS_FLOAT
#ifdef OMP
#pragma omp critical (lapack)
#endif
  ssyevr_(&JOBZ, &RANGE, &UPLO, &N,
			 A[0], &LDA, &VL, &VU,
			 &IL, &IU, &ABSTOL,
			 &M, D, Vt[0], &LDZ,
			 ISUPPZ, WORK, &LWORK,
			 IWORK, &LIWORK, &INFO);  // V is returned as V'!!!
#endif
#endif
  //  unsigned int rows = problemsize, columns = problemsize;
/*    unsigned int m = rows, n = columns, LDA = rows, LDU = rows, LDVT = columns; */
/*    unsigned int mindim = uimin(rows, columns); */
/*    unsigned int LWORK = uimax(3*uimin(m,n)+uimax(m,n),5*uimin(m,n)); */
/*    int INFO; */
/*    Matrix Xcolumnmajor = Matrix_allocate(columns, rows); */
/*    Matrix UT = Matrix_allocate(rows, rows); */
/*    Matrix temp = Matrix_allocate(columns, rows); */
/*    real D[mindim], WORK[LWORK], tol; */
/*    unsigned int i, j; */

/*    for (i = 0; i < rows; i++) */
/*       for (j = 0; j < columns; j++) */
/*          Xcolumnmajor[j][i] = .5 * (Lhat[i][j] + Lhat[j][i]); */
/* 	//switched to above from just Lhat[i][j] */

/* #ifdef REAL_IS_DOUBLE */
/* #ifdef OMP */
/* #pragma omp critical (lapack) */
/* #endif */
/* #ifdef _AIX */
/*    dgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D, */
/*            UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO); */
/* #else */
/*    dgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D, */
/*            UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO); */
/* #endif */
/* #else */
/* #ifdef REAL_IS_FLOAT */
/* #ifdef OMP */
/* #pragma omp critical (lapack) */
/* #endif */
/* #ifdef _AIX */
/*    sgesvd(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D, */
/*            UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO); */
/* #else */
/*    sgesvd_(&JOBU, &JOBVT, &m, &n, Xcolumnmajor[0], &LDA, D, */
/*            UT[0], &LDU, V[0], &LDVT, WORK, &LWORK, &INFO); */
/* #endif */
/* #endif */
/* #endif */

/*    if (INFO) { */
/*       printf("s/dgesvd returned error code %d\n", INFO); */
/*       exit(-2); */
/*    } */

/*    Matrix_free(Xcolumnmajor); */


	printf("trying to set penalty tol\n");
	if (penaltyType == ABSOLUTE_TOL) {
	  printf("absolute tol\n");
	  tol = tolerance;
	  *maxEig = D[0];
	} else if (penaltyType == RELATIVE_LHAT_TOL) {
	  printf("relative lhat tol\n");
	  tol = D[0] * tolerance;  // should be an option: abs or rel
	  *maxEig = D[0];
	} else if (penaltyType == RELATIVE_RAYLEIGH_TOL) {
	///////
	printf("relative rayleigh tol\n");
	  Vector  x = Vector_allocate(bound->numvariablecharges);
	  Vector Lx = Vector_allocate(bound->numvariablecharges);
	  for (i=0; i < bound->numvariablecharges; i++)
		 x[i] = Vt[N-1][i];
	  
	  Optimizer_multiplyByL(bound, unbound, Lx, x);
	  *maxEig =  Vector_dot(Lx, x, bound->numvariablecharges);
	  tol = tolerance * *maxEig;
	  printf("maxEig est = %f, tol = %f * maxEig = %f", *maxEig, tolerance, tol);
	  Vector_free(x);
	  Vector_free(Lx);
	}
	unsigned int dropCount = 0;
	for (i = 0; i < N; i++) {
	  if (D[i] <= tol) {
		 S[i][i] = (real) penaltyScale;
		 printf("penalizing vector %d\n", i);
		 dropCount++;
	  } else {
		 printf("%f > tol\n", D[i]);
		 S[i][i] = 0.0;
		  }
	}
	printf("penalizing %d-dimensional eigenspace...\n", dropCount);
	Matrix_copy(V, Vt, N, N);
	Matrix_transpose(&V, N, N);
	Matrix_copy(leftSingularVectors, V, N, N);
	Matrix_copy(rightSingularVectors, Vt, N, N);
	// see dgesvd man page.  V is actually returned VT... hence, we
	// DON'T transpose the returned variable
	
   Matrix_multiplymatrix(temp, V, S, N, N, N);  // no longer transposed
   Matrix_multiplymatrix(penalty, temp, Vt, N, N, N); // now "transposed" in Matrix structure

	for (i = 0; i < N; i++)
	  singularValues[i] = D[i];
	
   Matrix_free(V);
   Matrix_free(Vt);
   Matrix_free(A);
   Matrix_free(S);
	free(IWORK);
	free(ISUPPZ);
	free(WORK);
	Vector_free(D);
   Matrix_free(temp);
/* 	Matrix_writefile("penalty", penalty, N, N); */
/* 	Matrix_writefile("U", leftSingularVectors, N, N); */
/* 	Matrix_writefile("V", rightSingularVectors, N, N); */
}

void Optimizer_computePreconditionerHessian(PBEproblem bound, PBEproblem unbound, Matrix Lhat) {
   unsigned int i, j;
   Vector q = Vector_allocate(bound->numvariablecharges);
   Vector Lhatb_i = Vector_allocate(bound->numvariablecharges);
   Vector Lhatu_i = Vector_allocate(bound->numvariablecharges);
   Vector phi_b = Vector_allocate(bound->numtotalsurfacevariables);
   Vector phi_u = Vector_allocate(unbound->numtotalsurfacevariables);
   Vector Pphi_b = Vector_allocate(bound->numtotalsurfacevariables);
   Vector Pphi_u = Vector_allocate(unbound->numtotalsurfacevariables);

   for (i = 0; i < bound->numvariablecharges; i++) {
      Vector_zero(q, bound->numvariablecharges);
      q[i] = 1.0;

		PBEproblem_applyA1(bound, phi_b, q);
      PBEproblem_applyPreconditioner(bound, Pphi_b, phi_b);
      PBEproblem_applyA3(bound, Lhatb_i, Pphi_b);
		printf("phi_b[0] = %f  Pphi_b[0] = %f  Lhatb_i[0] = %f\n", phi_b[0], Pphi_b[0], Lhatb_i[0]);

      PBEproblem_applyA1(unbound, phi_u, q);
      PBEproblem_applyPreconditioner(unbound, Pphi_u, phi_u);
      PBEproblem_applyA3(unbound, Lhatu_i, Pphi_u);
		printf("phi_u[0] = %f  Pphi_u[0] = %f  Lhatu_i[0] = %f\n", phi_u[0], Pphi_u[0], Lhatu_i[0]);

      for (j = 0; j < bound->numvariablecharges; j++)
		  Lhat[j][i] = (Lhatb_i[j] - Lhatu_i[j]);
   }
	printf("about to write Lhat\n");
	//	Matrix_writefile("Lhat.m", Lhat, bound->numvariablecharges,  bound->numvariablecharges);
   Vector_free(q);
   Vector_free(phi_b);
   Vector_free(phi_u);
   Vector_free(Pphi_b);
   Vector_free(Pphi_u);
   Vector_free(Lhatb_i);
   Vector_free(Lhatu_i);
}

void Optimizer_updateChargeDistributions(PBEproblem bound, PBEproblem unbound) {
  PDBentry* PDBentries = bound->pdbentries;
  unsigned int numPDBentries = bound->numpdbentries;
  unsigned int i, counter = 0;
  Charge_free(bound->qualocationoperator->charges);
  Charge_free(unbound->qualocationoperator->charges);

  bound->qualocationoperator->charges = Charge_allocate();
  bound->qualocationoperator->charges->numcharges = numPDBentries;
  bound->qualocationoperator->charges->points = (Vector3D*)calloc(numPDBentries, sizeof(Vector3D));
  bound->qualocationoperator->charges->charges = Vector_allocate(numPDBentries);
  for (i = 0; i < numPDBentries; i++) {
	 Vector3D charge = Vector3D_allocate();
	 charge->x = PDBentries[i].x;
	 charge->y = PDBentries[i].y;
	 charge->z = PDBentries[i].z;
	 bound->qualocationoperator->charges->points[i] = charge;
	 bound->qualocationoperator->charges->charges[i] = PDBentries[i].charge;
  }
  
  unbound->qualocationoperator->charges = Charge_allocate();
  unbound->qualocationoperator->charges->points = (Vector3D*)calloc(numPDBentries, sizeof(Vector3D));
  unbound->qualocationoperator->charges->numcharges = unbound->numvariablecharges + unbound->numfixedligandcharges;
  unbound->qualocationoperator->charges->charges = Vector_allocate(numPDBentries);
  for (i = 0; i < numPDBentries; i++) {
	 if (toupper(PDBentries[i].chain) == 'V') {
		printf("%d matches to %d\n", i, counter);
		Vector3D charge = Vector3D_allocate();
		charge->x = PDBentries[i].x;
		charge->y = PDBentries[i].y;
		charge->z = PDBentries[i].z;
		unbound->qualocationoperator->charges->points[counter] = charge;
		unbound->qualocationoperator->charges->charges[counter] = PDBentries[i].charge;
		counter++;
	 }
  }
  
  if (usequalocation) {
	 Tree_free(bound->qualocationoperator->M1M3);
	 generateQualocationOperatorA1A3(bound->qualocationoperator);
	 Tree_free(unbound->qualocationoperator->M1M3);
	 generateQualocationOperatorA1A3(unbound->qualocationoperator);
  } else {
	 printf("Optimizer_updateChargeDistributions only works for qualocation!");
	 exit(-1);
  }
  
}

void Optimizer_computeLinearTerm(PBEproblem bound, PBEproblem unbound, Vector c) {
   unsigned int i;
   Vector boundPhiReact = Vector_allocate(bound->numvariablecharges);
   Vector unboundPhiReact = Vector_allocate(bound->numvariablecharges);
   Vector boundCoulombicTerm = Vector_allocate(bound->numvariablecharges);
   Vector unboundCoulombicTerm = Vector_allocate(bound->numvariablecharges);
	printf("in Optimizer_computeLinearTerm\n");
   Vector_zero(c, bound->numvariablecharges);
   
   Optimizer_zeroChainCharges(bound, bound->globalCharges);
   PBEproblem_solve(bound);
   PBEproblem_getVariableReactionPotentials(bound, boundPhiReact);
   Vector_addvector(c, boundPhiReact, bound->numvariablecharges);

   Optimizer_zeroChainCharges(unbound, unbound->globalCharges);
   PBEproblem_solve(unbound);
   PBEproblem_getVariableReactionPotentials(unbound, unboundPhiReact);
   Vector_addscaledvector(c, -1.0, unboundPhiReact, bound->numvariablecharges);

/* 	Vector_writefile("opt_boundSol", bound->Sol, bound->numtotalsurfacevariables); */
/* 	Vector_writefile("opt_unboundSol", unbound->Sol, unbound->numtotalsurfacevariables); */
/* 	Vector_writefile("opt_boundRHS", bound->RHS, bound->numtotalsurfacevariables); */
/* 	Vector_writefile("opt_unboundRHS", unbound->RHS, unbound->numtotalsurfacevariables); */

/*    for (i = 0; i < unbound->numvariablecharges; i++) */
/*       printf("reactphidiff[i] = %f\n", c[i]); */
   
   Optimizer_getCoulombicInteraction(bound, boundCoulombicTerm);
   Vector_addvector(c, boundCoulombicTerm, bound->numvariablecharges);
   Optimizer_getCoulombicInteraction(unbound, unboundCoulombicTerm);
   Vector_addscaledvector(c, -1.0, unboundCoulombicTerm, bound->numvariablecharges);
    
   Vector_free(boundCoulombicTerm);
   Vector_free(unboundCoulombicTerm);
   Vector_free(boundPhiReact);
   Vector_free(unboundPhiReact);
	printf("exiting Optimizer_computeLinearTerm\n");
}

void Optimizer_zeroChainCharges(PBEproblem problem, Vector fixedCharges) {
   unsigned int i;
   
   for (i = 0; i < problem->numpdbentries; i++) {
      if (variablechain == problem->pdbentries[i].chain)
         fixedCharges[i] = 0.0;
      else
         fixedCharges[i] = problem->pdbentries[i].charge;
   }
}

void Optimizer_getCoulombicInteraction(PBEproblem problem, Vector interaction) {
   unsigned int i, j;
   Vector3D *location;
   location = (Vector3D *)calloc(problem->numpdbentries, sizeof(Vector3D));
   
   Vector_zero(problem->globalPhiReact, problem->numpdbentries);

   for (i = 0; i < problem->numpdbentries; i++) {
      location[i] = Vector3D_allocate();
      location[i]->x = problem->pdbentries[i].x;
      location[i]->y = problem->pdbentries[i].y;
      location[i]->z = problem->pdbentries[i].z;
   }

   for (i = 0; i < problem->numpdbentries; i++) {
      for (j = 0; j < problem->numpdbentries; j++) {
         if (i == j)
            continue;
         problem->globalPhiReact[i] += .592 * KT_CONVERSION * problem->pdbentries[j].charge /
            (Vector3D_distance(location[i], location[j]) * innerdielectric);
      }
   }

   for (i = 0; i < problem->numpdbentries; i++)
      Vector3D_free(location[i]);
   free(location);
   
   PBEproblem_getVariableReactionPotentials(problem, interaction);
}

void Optimizer_multiplyByL(PBEproblem bound, PBEproblem unbound, Vector Lx, Vector x) {
  unsigned int n_c = bound->numvariablecharges;
  Vector unboundReactPot = Vector_allocate(n_c);
  Vector boundReactPot = Vector_allocate(n_c);
  PBEproblem_setVariableChargeVector(unbound, x);
  printf("about to solve unbound\n");
  PBEproblem_solve(unbound);
  PBEproblem_getVariableReactionPotentials(unbound, unboundReactPot);

  num_GMRES_iter = 0;
  PBEproblem_setVariableChargeVector(bound, x);
  printf("about to solve bound\n");
  PBEproblem_solve(bound);
  PBEproblem_getVariableReactionPotentials(bound, boundReactPot);
  
  Vector_copy(Lx, boundReactPot, n_c);
  Vector_subtractvector(Lx, unboundReactPot, n_c);

  Vector_free(unboundReactPot);
  Vector_free(boundReactPot);
}
