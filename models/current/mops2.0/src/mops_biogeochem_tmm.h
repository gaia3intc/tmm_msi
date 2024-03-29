/* $Header: /Users/ikriest/CVS/mops/mops_biogeochem.h,v 1.2 2015/11/17 14:18:51 ikriest Exp $ */
/* $Name: mops-2_0 $ */
/* T.Tanioka added fbgc8, fbgc9, fbgc10, fbgc 11 (Nov 2020) */
/* T.Tanioka added READ_MARTINB (Dec 2020) */
/* K.Matsumoto added fbgc 12 (PP of PFT2) and 13 (C:P of PFT2) (Aug, 2021) */

extern void mops_biogeochem_copy_data_(PetscInt *nzloc, PetscInt *itr, PetscScalar localTR[], PetscScalar localJTR[], 
                                PetscScalar *DeltaT, PetscInt *direction);

extern void insolation_(PetscInt *N, PetscScalar *myTime, PetscScalar locallatitude[], PetscScalar *daysperyear, PetscScalar localswrad[], PetscScalar localstau[]);

extern void mops_biogeochem_ini_(PetscInt *Nrloc, PetscScalar *DeltaT, 
#ifdef CARBON                      
                                   PetscScalar *localph,
#endif
                                   PetscScalar localTs[], PetscScalar localSs[], 
                                   PetscScalar localdz[], PetscScalar drF[], PetscInt *nzmax, PetscInt *nzeuph,
                                   PetscInt *numBiogeochemStepsPerOceanStep,
                                   PetscBool *setDefaults);

extern void mops_biogeochem_model_(PetscInt *Nrloc, PetscScalar *DeltaT,
#ifdef CARBON                      
				   PetscScalar *DICglobavg, PetscScalar *ALKglobavg, PetscScalar *localEmP, PetscScalar *localpCO2atm,
#endif
                                   PetscScalar localTs[],PetscScalar localSs[], 
                                   PetscScalar *localfice, PetscScalar *localswrad, PetscScalar *localstau,
                                   PetscScalar *localwind, PetscScalar *localatmosp, PetscScalar localdz[], 
#ifdef READ_MARTINB
                                   PetscScalar *localmartinbc, PetscScalar drF[], PetscInt *nzmax,
#endif
#ifdef CARBON                      
                                   PetscScalar *localph,
				   PetscScalar *localco2airseaflux,
#endif
                                   PetscScalar *localburial, PetscScalar *GRunoff, PetscScalar localrunoffvol[],
                                   PetscBool *useSeparateBiogeochemTimeStepping);
/*#ifndef ORGCARBON
/*extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
/*                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
/*					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
/*                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[]);
/*#else
/*
/*#ifndef PFT
/* fbgc9 = Sediment_C, fbgc10 = Phytoplankton C:P uptake ratio, fbgc11 = Zooplankton C:P uptake ratio */ 
/*extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
/*                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
/*					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
/*                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[],
/*                                         PetscScalar localfbgc10[], PetscScalar localfbgc11[]);
/*#else
/* fbgc12 = PP and fbgc13 = C:P uptake ratio of second PFT */
/*extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
/*                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
/*					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
/*                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[],
/*                                         PetscScalar localfbgc10[], PetscScalar localfbgc11[],
/*                                         PetscScalar localfbgc12[], PetscScalar localfbgc13[]);
/*#endif
/*
/*#endif*/

#ifndef PFT
#ifdef ORGCARBON
#ifndef FLEXCP
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[]);
#else
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[],
                                         PetscScalar localfbgc10[], PetscScalar localfbgc11[]);
#endif

#else
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[]);
#endif 
      
#else

#ifndef ORGCARBON      
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[]);
#else

#ifndef FLEXCP      
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[],
                                         PetscScalar localfbgc10[]);
#else
extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], 
                                         PetscScalar localfbgc7[], PetscScalar localfbgc8[], PetscScalar localfbgc9[],
                                         PetscScalar localfbgc10[], PetscScalar localfbgc11[],
                                         PetscScalar localfbgc12[], PetscScalar localfbgc13[]);
#endif
#endif      
#endif      

extern void mops_biogeochem_set_params_(PetscInt *numbgcparams, PetscScalar bgcparams[]);

extern void mops_biogeochem_misfit_(PetscInt *Nrloc,
                                         PetscScalar localmbgc1[], PetscScalar localmbgc2[], PetscScalar localmbgc3[]);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define mops_biogeochem_copy_data_ mops_biogeochem_copy_data
#define insolation_ insolation
#define mops_biogeochem_ini_ mops_biogeochem_ini
#define mops_biogeochem_model_ mops_biogeochem_model
#define mops_biogeochem_diagnostics_ mops_biogeochem_diagnostics
#define mops_biogeochem_set_params_ mops_biogeochem_set_params
#define mops_biogeochem_misfit_ mops_biogeochem_misfit
#endif 

