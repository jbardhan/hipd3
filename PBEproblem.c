#include "PBEproblem.h"
#include "Overlap.h"
#include "Integration.h"
#include <sys/time.h>
#include <sys/resource.h>

PBEproblem PBEproblem_allocate(char *PDBfilename, char *SRFfilename, int compress) {
  PBEproblem problem = NULL;
  problem = (PBEproblem)(calloc(1, sizeof(_PBEproblem)));
  PBEproblem_allocateButDontCompress(problem, PDBfilename, SRFfilename);
  if (compress)
	 PBEproblem_initialize(problem);

  return problem;
}

PBEproblem PBEproblem_loadOnlyPDB(char *PDBfilename) {
  PBEproblem problem = NULL;
  problem = (PBEproblem)(calloc(1, sizeof(_PBEproblem)));
  problem->preconditioner = NULL;
  problem->qualocationoperator = NULL;
  problem->pbesurfaceoperator = NULL;
  readPDB(PDBfilename, &(problem->numpdbentries), &(problem->pdbentries));
  assignRadiiCharges(problem->pdbentries, problem->numpdbentries, SIZentries, numSIZentries,
                      CRGentries, numCRGentries);

  return problem;
}

void PBEproblem_allocateButDontCompress(PBEproblem problem, char *PDBfilename, char *SRFfilename) {
	problem->preconditioner = NULL;
	problem->qualocationoperator = NULL;
	problem->pbesurfaceoperator = NULL;
   problem->useOverlap = 0; // default preconditioner is the diagonal
									 // or block-diagonal; set manually to get
									 // overlap preconditioner!
   readPDB(PDBfilename, &(problem->numpdbentries), &(problem->pdbentries));
   assignRadiiCharges(problem->pdbentries, problem->numpdbentries, SIZentries, numSIZentries,
                      CRGentries, numCRGentries);

   readSRF(SRFfilename,
           &(problem->saltpanels), &(problem->numsaltpanels), &(problem->numsalts),
           &(problem->dielectricpanels), &(problem->numdielectricpanels),
           &(problem->numdielectrics),  &(problem->dielectricparent),
           &(problem->dielectriccavitypanels), &(problem->numdielectriccavitypanels),
           &(problem->numdielectriccavities), &(problem->dielectriccavityparent),
           &(problem->saltcavitypanels), &(problem->numsaltcavitypanels),
           &(problem->numsaltcavities), &(problem->saltcavityparent),
           &(problem->numtotalpanels));

   problem->numvariablecharges = 0;
   PBEproblem_findVariableChargeIndices(problem, variablechain);
   problem->numfixedligandcharges = 0;
   PBEproblem_findFixedLigandChargeIndices(problem, fixedligandchain);
   problem->numfixedreceptorcharges = 0;
   PBEproblem_findFixedReceptorChargeIndices(problem, variablechain, fixedligandchain);
   
   problem->globalCharges = Vector_allocate(problem->numpdbentries);

   if (usequalocation)
      problem->numtotalsurfacevariables = problem->numtotalpanels;
   else
      problem->numtotalsurfacevariables = 2 * problem->numtotalpanels;

   problem->RHS = Vector_allocate(problem->numtotalsurfacevariables);
   problem->Sol = Vector_allocate(problem->numtotalsurfacevariables);

   problem->globalPhiReact = Vector_allocate(problem->numpdbentries);
   problem->variablePhiReact = Vector_allocate(problem->numvariablecharges);
}

void PBEproblem_initialize(PBEproblem problem)
{
   printf("found %d %d %d charges\n",
          problem->numvariablecharges,
          problem->numfixedligandcharges,
          problem->numfixedreceptorcharges);

   struct rusage ruse;
   struct timeval tval;
   real starttime, endtime, startwalltime, endwalltime;

   if (usequalocation) {
      printf("generating qualocation operator\n");
   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
   generateQualocationOperator(&(problem->qualocationoperator),
                               problem->pdbentries, problem->numpdbentries,
                               problem->dielectricpanels, problem->numdielectricpanels, problem->numdielectrics,
                               problem->dielectriccavitypanels, problem->numdielectriccavitypanels, problem->numdielectriccavities,
                               problem->dielectriccavityparent, problem->numtotalpanels);
   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;
   printf("time to initialize operator: %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);

   generateQualocationOperatorPreconditioner(&problem->preconditioner,
                                             problem->qualocationoperator,
                                             problem->numtotalpanels);
   } else {
      printf("generating Green's theorem operator\n");
      
   generateSurfaceOperator(&(problem->pbesurfaceoperator), problem->pdbentries, problem->numpdbentries,
                           problem->saltpanels, problem->numsaltpanels, problem->numsalts,
                           problem->dielectricpanels, problem->numdielectricpanels,
                           problem->numdielectrics, problem->dielectricparent,
                           problem->dielectriccavitypanels, problem->numdielectriccavitypanels,
                           problem->numdielectriccavities, problem->dielectriccavityparent,
                           problem->saltcavitypanels, problem->numsaltcavitypanels,
                           problem->numsaltcavities, problem->saltcavityparent,
                           problem->numtotalpanels);
   generateSurfaceOperatorPreconditioner(&(problem->preconditioner), problem->pbesurfaceoperator, problem->numtotalpanels);
   }
}


void PBEproblem_generateLhatOnTheFly(PBEproblem problem, Matrix Lhat) {
/*   generateQualocationOperatorMinimal(&(problem->qualocationoperator), */
/*                                problem->pdbentries, problem->numpdbentries, */
/*                                problem->dielectricpanels, problem->numdielectricpanels, problem->numdielectrics, */
/*                                problem->dielectriccavitypanels, problem->numdielectriccavitypanels, problem->numdielectriccavities, */
/*                                problem->dielectriccavityparent, problem->numtotalpanels); */

/*   problem->preconditioner = NULL; */
/*   PBEproblem_generateDiagonalPreconditionerOnTheFly(problem, &(problem->preconditioner)); */
  printf("generateLhatOnTheFly is commented out.\n");exit(-1);
}

void PBEproblem_generateDiagonalPreconditionerOnTheFly(PBEproblem problem, Preconditioner *P) {
  unsigned int i, j, count = 0;
  if (*P != NULL) {
	 Preconditioner_free(*P);
  }
  *P = Preconditioner_allocate(problem->numtotalsurfacevariables, problem->numtotalsurfacevariables);

  if (usequalocation) { // basically Preconditioner_fill_diagonal_solv_ecf_qual without Tree_extractdiagonal
	 Vector diag = Vector_allocate(problem->numtotalpanels);

	 Panel* allpanels = (Panel *)calloc(problem->numtotalpanels, sizeof(Panel));
	 
	 for (i = 0; i < problem->numdielectrics; i++)
      for (j = 0; j < problem->numdielectricpanels[i]; j++) {
         allpanels[count] = problem->dielectricpanels[i][j];
         count++;
      }

   for (i = 0; i < problem->numdielectriccavities; i++) {
      if (problem->dielectriccavityparent[i] == 999)
         continue;

      for (j = 0; j < problem->numdielectriccavitypanels[i]; j++) {
         allpanels[count] = problem->dielectriccavitypanels[i][j];
         count++;
      }
   }
	real idiel = innerdielectric;
	real odiel = outerdielectric;
	for (i = 0; i < problem->numtotalpanels; i++) {
	  diag[i] = Integration(allpanels[i]->centroid, allpanels[i], POISSON_KERNEL, NULL, DOUBLE_LAYER_INT);
	  Preconditioner_set(*P, i, i, (-odiel / ((odiel - idiel) * idiel) + (diag[i]) / (4.0 * M_PI * idiel)) * allpanels[i]->area);
	}

	Vector_free(diag);
  } else {
	 printf("PBEproblem_generateDiagonalPreconditionerOnTheFly: Green's theorem case not implemented yet!\n");
  }
  Preconditioner_factor(*P);
}

void PBEproblem_loadNewChargeDistribution(PBEproblem problem, char *PDBfilename, char *CRGfilename) {
   readCRG(CRGfilename, &numCRGentries, &CRGentries);

	printf("loaded CRG\n");
   readPDB(PDBfilename, &(problem->numpdbentries), &(problem->pdbentries));
	printf("loaded new PDB\n");
   assignRadiiCharges(problem->pdbentries, problem->numpdbentries, SIZentries, numSIZentries,
                      CRGentries, numCRGentries);
	printf("assignRadiiCharges has returned\n");
   problem->numvariablecharges = 0;
   PBEproblem_findVariableChargeIndices(problem, variablechain);
   problem->numfixedligandcharges = 0;
   PBEproblem_findFixedLigandChargeIndices(problem, fixedligandchain);
   problem->numfixedreceptorcharges = 0;
   PBEproblem_findFixedReceptorChargeIndices(problem, variablechain, fixedligandchain);
	printf("indices updated\n");

	free(problem->globalCharges);
	free(problem->globalPhiReact);
	free(problem->variablePhiReact);

	problem->globalCharges = Vector_allocate(problem->numpdbentries);
	problem->globalPhiReact = Vector_allocate(problem->numpdbentries);
	problem->variablePhiReact = Vector_allocate(problem->numvariablecharges);

}

void PBEproblem_findVariableChargeIndices(PBEproblem problem, char variablechain) {
   unsigned int i, currentcount = 0;
   problem->numvariablecharges = 0;
   for (i = 0; i < problem->numpdbentries; i++)
      if (variablechain == problem->pdbentries[i].chain)
         problem->numvariablecharges++;

   problem->variablechargeindextoglobalindex = (unsigned int *)calloc(problem->numvariablecharges, sizeof(unsigned int));

   for (i = 0; i < problem->numpdbentries; i++) {
      if (variablechain == problem->pdbentries[i].chain) {
         problem->variablechargeindextoglobalindex[currentcount] = i;
         currentcount++;
      }
   }
}

void PBEproblem_findFixedLigandChargeIndices(PBEproblem problem, char fixedligandchain) {
   unsigned int i, currentcount = 0;
   for (i = 0; i < problem->numpdbentries; i++)
      if (fixedligandchain == problem->pdbentries[i].chain)
         problem->numfixedligandcharges++;

   problem->fixedligandchargeindextoglobalindex = (unsigned int *)calloc(problem->numfixedligandcharges, sizeof(unsigned int));

   for (i = 0; i < problem->numpdbentries; i++) {
      if (fixedligandchain == problem->pdbentries[i].chain) {
         problem->fixedligandchargeindextoglobalindex[currentcount] = i;
         currentcount++;
      }
   }
}

void PBEproblem_findFixedReceptorChargeIndices(PBEproblem problem, char variablechain, char fixedligandchain) {
   unsigned int i, currentcount = 0;
   for (i = 0; i < problem->numpdbentries; i++)
	  if (! ((variablechain == problem->pdbentries[i].chain) || (fixedligandchain==problem->pdbentries[i].chain)) )
         problem->numfixedreceptorcharges++;

   problem->fixedreceptorchargeindextoglobalindex = (unsigned int *)calloc(problem->numfixedreceptorcharges, sizeof(unsigned int));

   for (i = 0; i < problem->numpdbentries; i++) {
	  if (! ((variablechain == problem->pdbentries[i].chain) || (fixedligandchain==problem->pdbentries[i].chain)) ) {
         problem->fixedreceptorchargeindextoglobalindex[currentcount] = i;
         currentcount++;
      }
   }
}

void PBEproblem_free(PBEproblem problem) {
   if (usequalocation) {
	  if (problem->qualocationoperator)
      QualocationOperator_free(problem->qualocationoperator);
   } else {
	  if (problem->pbesurfaceoperator)
      SurfaceOperator_free(problem->pbesurfaceoperator);
   }
	if (problem->preconditioner)
	  Preconditioner_free(problem->preconditioner);
	
   // what about freeing panels, PDBentries?
   free(problem->variablechargeindextoglobalindex);

   Vector_free(problem->globalCharges);
   Vector_free(problem->RHS);
   Vector_free(problem->Sol);
   Vector_free(problem->globalPhiReact);
   Vector_free(problem->variablePhiReact);
   free(problem);
}

void PBEproblem_setVariableChargeVector(PBEproblem problem, Vector variablechargevec) {
   unsigned int i;
	if (problem->globalCharges == variablechargevec) {
	  return;
	}
   // assumes arguments are appropriately dimensioned
   for (i = 0; i < problem->numpdbentries; i++)
      problem->globalCharges[i] = 0.0;
   
   for (i = 0; i < problem->numvariablecharges; i++) {
      problem->globalCharges[problem->variablechargeindextoglobalindex[i]] = variablechargevec[i];
   }
}

void PBEproblem_getVariableReactionPotentials(PBEproblem problem, Vector localreactionvector) {
   unsigned int i;
   for (i = 0; i < problem->numvariablecharges; i++) 
      localreactionvector[i] = problem->globalPhiReact[problem->variablechargeindextoglobalindex[i]];
}

void PBEproblem_solve(PBEproblem problem) {
   unsigned int i;
   struct rusage ruse;
   struct timeval tval;
   real starttime, endtime, startwalltime, endwalltime;
   getrusage(RUSAGE_SELF, &ruse);
   starttime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   startwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

   Vector_zero(problem->RHS, problem->numtotalsurfacevariables);
   Vector_zero(problem->Sol, problem->numtotalsurfacevariables);
   if (usequalocation) {
      QualocationOperator_makeRHS_fromVector(problem->qualocationoperator,
                                             problem->RHS, problem->globalCharges, problem->numpdbentries);
      //			Vector_writefile("b_qual_known.m", problem->RHS, problem->numtotalsurfacevariables);

/* 		printf("pbeproblem_solve_RHS = ["); */
/* 		for (i = 0; i < problem->numtotalsurfacevariables; i++) { */
/* 		  printf("%f\n", problem->RHS[i]); */
/* 		} */
/* 		printf("];\n"); */

      GMRES_QualocationOperator(problem->qualocationoperator, problem->preconditioner,
                                problem->RHS, problem->Sol,
                                problem->numtotalsurfacevariables, errortolerance);

      //		Vector_writefile("x_qual_known.m",problem->Sol,problem->numtotalsurfacevariables);
/* 		printf("pbeproblem_solve_Sol = ["); */
/* 		for (i = 0; i < problem->numtotalsurfacevariables; i++) { */
/* 		  printf("%f\n", problem->Sol[i]); */
/* 		} */
/* 		printf("];\n"); */

QualocationOperator_collectPotentials_toVector(problem->qualocationoperator,
                                                     problem->globalPhiReact, problem->Sol);

   } else {
      SurfaceOperator_makeRHS_fromVector(problem->pbesurfaceoperator, problem->RHS,
                                         problem->globalCharges, problem->numpdbentries);
      //		Vector_writefile("b_yoon_coll.m", problem->RHS, problem->numtotalsurfacevariables);
/* 		printf("rhs=["); */
/* 		for (i = 0; i < problem->numtotalsurfacevariables; i++) { */
/* 		  printf("%f\n", problem->RHS[i]); */
/* 		} */
/* 		printf("];\n"); */
      GMRES_SurfaceOperator(problem->pbesurfaceoperator, problem->preconditioner,
                            problem->RHS, problem->Sol,
                            problem->numtotalsurfacevariables, errortolerance);
      Vector_zero(problem->globalPhiReact, problem->numpdbentries);
      SurfaceOperator_collectPotentials_toVector(problem->pbesurfaceoperator,
                                                 problem->globalPhiReact, problem->Sol);
      //		Vector_writefile("x_yoon_coll.m", problem->Sol, problem->numtotalsurfacevariables);
/* 		printf("x=["); */
/* 		for (i = 0; i < problem->numtotalsurfacevariables; i++) { */
/* 		  printf("%f\n", problem->Sol[i]); */
/* 		} */
/* 		printf("];\n"); */
   }
   getrusage(RUSAGE_SELF, &ruse);
   endtime = ruse.ru_utime.tv_sec + ruse.ru_stime.tv_sec +
      1e-6 * (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec);
   gettimeofday(&tval, NULL);
   endwalltime = tval.tv_sec + 1e-6 * tval.tv_usec;

   printf("time for PBEproblem_solve(): %.2f s (%.2f)\n", endtime - starttime, endwalltime - startwalltime);

   for (i = 0; i < problem->numpdbentries; i++)
      problem->globalPhiReact[i] *= .592 * KT_CONVERSION / innerdielectric;
}

void PBEproblem_applyPreconditioner(PBEproblem problem, Vector Px, Vector x) {
  if (problem->useOverlap) {
	 Preconditioner_multiply(Px, problem->overlapPreconditioner, x, problem->numtotalsurfacevariables);
  } else {
   Preconditioner_solve(Px, problem->preconditioner, x);
  }
}

void PBEproblem_applyA1(PBEproblem problem, Vector A1x, Vector x) {

  PBEproblem_setVariableChargeVector(problem, x);

  if (usequalocation) {
	 QualocationOperator_makeRHS_fromVector(problem->qualocationoperator,
														 A1x, problem->globalCharges,
														 problem->numpdbentries);
  } else {
	 SurfaceOperator_makeRHS_fromVector(problem->pbesurfaceoperator,
													A1x, problem->globalCharges,
													problem->numpdbentries);
  }
  // need to negate the calculated vector?

}

void PBEproblem_applyA2(PBEproblem problem, Vector A2x, Vector x) {
   if (usequalocation) {
      QualocationOperator_multiply(problem->qualocationoperator, A2x, x);
   } else {
      SurfaceOperator_topmultiply(problem->pbesurfaceoperator, A2x, x);
   }
}

void PBEproblem_applyA3(PBEproblem problem, Vector A3x, Vector x) {
   unsigned int i;
   Vector_zero(problem->globalPhiReact, problem->numpdbentries);

   if (usequalocation) {
      QualocationOperator_collectPotentials_toVector(problem->qualocationoperator,
                                                     problem->globalPhiReact, x);
   } else {
      SurfaceOperator_collectPotentials_toVector(problem->pbesurfaceoperator,
                                                 problem->globalPhiReact, x);
   }
   
   for (i = 0; i < problem->numpdbentries; i++)
      problem->globalPhiReact[i] *= .592 * KT_CONVERSION / innerdielectric;

   PBEproblem_getVariableReactionPotentials(problem, A3x);
}

void PBEproblem_setupOverlapPreconditioner(PBEproblem problem) {
  problem->useOverlap = 1;
  problem->overlapPreconditioner = Preconditioner_allocate(problem->numtotalsurfacevariables,
																			  problem->numtotalsurfacevariables);
  Preconditioner_fill_overlap_ecf_qual_cav(problem->overlapPreconditioner, problem->qualocationoperator->tree,
														 problem->qualocationoperator->tree->panels,
														 problem->numtotalsurfacevariables, problem->numtotalsurfacevariables,
														 problem->qualocationoperator->innerdielectric,
														 problem->qualocationoperator->outerdielectric);
}


void PBEproblem_saveOverlapPreconditioner(PBEproblem problem, char *filename) {
  unsigned int i,j;
  if (problem->overlapPreconditioner == NULL) {
	 printf("ERROR: no overlap preconditioner has been created yet.\n");
	 return;
  }
  Vector Pqq = Vector_allocate(problem->numtotalsurfacevariables);
  Vector qq = Vector_allocate(problem->numtotalsurfacevariables);
  FILE *OUTFILE = fopen(filename, "w");
  for (i = 0; i < problem->numtotalsurfacevariables; i++) {
	 Vector_zero(qq, problem->numtotalsurfacevariables); qq[i]=1.0;
	 Preconditioner_multiply(Pqq, problem->overlapPreconditioner, qq, problem->numtotalsurfacevariables);
	 for (j = 0; j < problem->numtotalsurfacevariables; j++) {
		 fprintf(OUTFILE, "%f  ", Pqq[j]);
	 }
	 fprintf(OUTFILE, "\n");
  }
  fclose(OUTFILE);
  Vector_free(qq);
  Vector_free(Pqq);
}
