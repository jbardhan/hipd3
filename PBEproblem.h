#ifndef __PBEPROBLEM_H__
#define __PBEPROBLEM_H__

#include "FFTSVDpbeAPI.h"

typedef struct _PBEproblem {
   unsigned int numsalts;
   Panel** saltpanels;
   unsigned int* numsaltpanels;

   unsigned int numdielectrics;
   Panel** dielectricpanels;
   unsigned int* numdielectricpanels;
   unsigned int* dielectricparent; 

   unsigned int numdielectriccavities;
   Panel** dielectriccavitypanels;
   unsigned int* numdielectriccavitypanels;
   unsigned int* dielectriccavityparent;

   unsigned int numsaltcavities;
   Panel** saltcavitypanels;
   unsigned int* numsaltcavitypanels;
   unsigned int* saltcavityparent;

   unsigned int numtotalpanels;
   unsigned int numtotalsurfacevariables; //  == 2 * numtotalpanels in green's thm,
   // but only 1 * numtotalpanels in qualocation
   
   /* Surface Operator Related */
   SurfaceOperator pbesurfaceoperator;  /* The big momma */
   QualocationOperator qualocationoperator;  // either this or the above get used, not both!
   Preconditioner preconditioner;
   Preconditioner overlapPreconditioner;
   unsigned int useOverlap; /* for the moment, only implemented for
									   qualocation.  however, it has been
									   demonstrated to reduce the # of gmres
									   iterations by a factor of 3 for Yoon and
									   Lenhoff formulation*/

   PDBentry* pdbentries;
   unsigned int numpdbentries;

   unsigned int numvariablecharges;
   unsigned int *variablechargeindextoglobalindex;
   unsigned int numfixedligandcharges;
   unsigned int *fixedligandchargeindextoglobalindex;
   unsigned int numfixedreceptorcharges;
   unsigned int *fixedreceptorchargeindextoglobalindex;

   Vector globalCharges;
   Vector RHS;
   Vector Sol;
   Vector globalPhiReact;
   Vector variablePhiReact;

   char *directory;  /* the bound and unbound PDB files are in
                        separate directories -- each uses the
                        appropriate delphi.prl fort.10, etc files */
} _PBEproblem;

typedef struct _PBEproblem* PBEproblem;

PBEproblem PBEproblem_allocate(char *PDBfilename, char *SRFfilename, int compress);
PBEproblem PBEproblem_loadOnlyPDB(char *PDBfilename);
void PBEproblem_allocateButDontCompress(PBEproblem problem, char *PDBfilename, char *SRFfilename);
void PBEproblem_initialize(PBEproblem problem);
void PBEproblem_generateLhatOnTheFly(PBEproblem problem, Matrix Lhat);
void PBEproblem_generateDiagonalPreconditionerOnTheFly(PBEproblem problem, Preconditioner *P);
void PBEproblem_loadNewChargeDistribution(PBEproblem problem, char *PDBfilename, char *CRGfilename);
int PBEproblem_load(PBEproblem problem);
void PBEproblem_free(PBEproblem problem);
void PBEproblem_setVariableChargeVector(PBEproblem problem, Vector variableCharges);
void PBEproblem_getVariableReactionPotentials(PBEproblem problem, Vector variableReactPot);
void PBEproblem_solve(PBEproblem problem);
void PBEproblem_applyPreconditioner(PBEproblem problem, Vector Px, Vector x);
void PBEproblem_applyA1(PBEproblem problem, Vector A1x, Vector x);
void PBEproblem_applyA2(PBEproblem problem, Vector A2x, Vector x);
void PBEproblem_applyA3(PBEproblem problem, Vector A3x, Vector x);
void PBEproblem_findVariableChargeIndices(PBEproblem problem, char variablechain);
void PBEproblem_findFixedLigandChargeIndices(PBEproblem problem, char fixedligandchain);
void PBEproblem_findFixedReceptorChargeIndices(PBEproblem problem, char variablechain, char fixedligandchain);
void PBEproblem_saveOverlapPreconditioner(PBEproblem problem, char *filename);
void PBEproblem_setupOverlapPreconditioner(PBEproblem problem);

extern SIZentry* SIZentries;
extern unsigned int numSIZentries;
extern CRGentry* CRGentries;
extern unsigned int numCRGentries;
extern char variablechain;
extern char fixedligandchain;
extern char fixedreceptorchain;
   
#endif
