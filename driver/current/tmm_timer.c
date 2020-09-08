#define MAX_NUM_INTERVALS 1000

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_timer.h"

PetscErrorCode iniPeriodicTimer( const char pre[], PeriodicTimer *thetimer )
{

    PetscErrorCode ierr;
    PetscBool flg, flg1;
    PetscInt it;
	PetscViewer fd;
	int fp;
    char timeFile[PETSC_MAX_PATH_LEN];

/*  read time data */
    ierr = PetscOptionsGetReal(NULL,pre,"-cycle_period",&thetimer->cyclePeriod,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate cycling time with the -%scycle_period option",pre);

    ierr = PetscOptionsGetReal(NULL,pre,"-cycle_step",&thetimer->cycleStep,&flg);CHKERRQ(ierr);
    if (flg) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"WARNING!: Cycling step has been specified for periodic object %s\n",pre);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"  This is a legacy option retained for backward compatibility and will be removed in future releases\n");CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Use -%snum_per_period instead\n",pre);CHKERRQ(ierr);
      thetimer->numPerPeriod=thetimer->cyclePeriod/thetimer->cycleStep;	
/*    array for holding extended time array */
	  PetscMalloc((thetimer->numPerPeriod+2)*sizeof(PetscScalar), &thetimer->tdp);
	  for (it=0; it<=thetimer->numPerPeriod+1; it++) {
		thetimer->tdp[it]=(-thetimer->cycleStep/2.0) + it*thetimer->cycleStep;
	  } 	
    } else {
      ierr = PetscOptionsGetInt(NULL,pre,"-num_per_period",&thetimer->numPerPeriod,&flg1);CHKERRQ(ierr);
      if (!flg1) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of fields per period with the -%snum_per_period option",pre);
/*    array for holding extended time array */
	  PetscMalloc((thetimer->numPerPeriod+2)*sizeof(PetscScalar), &thetimer->tdp);
      ierr = PetscOptionsGetString(NULL,pre,"-periodic_times_file",timeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
      if (flg1) {
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,timeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		ierr = PetscBinaryRead(fp,&thetimer->tdp[1],thetimer->numPerPeriod,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		thetimer->tdp[0]=thetimer->tdp[thetimer->numPerPeriod]-thetimer->cyclePeriod;
		thetimer->tdp[thetimer->numPerPeriod+1]=thetimer->tdp[1]+thetimer->cyclePeriod;
      } else {
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"Assuming equally-spaced fields for periodic object %s\n",pre);CHKERRQ(ierr);            
        thetimer->cycleStep=thetimer->cyclePeriod/thetimer->numPerPeriod;
		for (it=0; it<=thetimer->numPerPeriod+1; it++) {
		  thetimer->tdp[it]=(-thetimer->cycleStep/2.0) + it*thetimer->cycleStep;
		} 	
      }
    }

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic object %s specified at times:\n",pre);CHKERRQ(ierr);            
	for (it=0; it<=thetimer->numPerPeriod+1; it++) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"tdp=%10.5f\n", thetimer->tdp[it]);CHKERRQ(ierr);        
	}    
    
    return 0;
    
}

PetscErrorCode iniStepTimer( const char pre[], PetscInt Iter0, StepTimer *thetimer )
{

    PetscErrorCode ierr;
    PetscBool flg;
    PetscInt it;
    PetscInt tmparr[MAX_NUM_INTERVALS];

/*  read time step data */
    thetimer->startTimeStep = Iter0 + 1; /* by default we start at first time step */
    ierr = PetscOptionsGetInt(NULL,pre,"-start_time_step",&thetimer->startTimeStep,&flg);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time step for StepTimer object %s is %d\n", pre, thetimer->startTimeStep);CHKERRQ(ierr);	  

    thetimer->maxNumIntervals=MAX_NUM_INTERVALS;
    ierr = PetscOptionsGetIntArray(NULL,pre,"-time_steps",tmparr,&thetimer->maxNumIntervals,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of step timer time steps with the -%stime_steps flag",pre);

    if (thetimer->maxNumIntervals==1) {
      thetimer->fixedStep=PETSC_TRUE;
	  thetimer->currInterval=0; /* Not used but we set it anyway to be safe */
      thetimer->numTimeSteps=tmparr[0];
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Fixed interval of %d specified for StepTimer object %s\n", thetimer->numTimeSteps, pre);CHKERRQ(ierr);	  
    } else {
      thetimer->fixedStep=PETSC_FALSE;      
	  PetscMalloc(thetimer->maxNumIntervals*sizeof(PetscInt), &thetimer->timeIntervals);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Variable number of intervals specified for StepTimer object %s\n", pre);CHKERRQ(ierr);	  
	  for (it=0; it<thetimer->maxNumIntervals; it++) {
		thetimer->timeIntervals[it] = tmparr[it];
		ierr = PetscPrintf(PETSC_COMM_WORLD,"  Interval #%d=%d\n", it+1,thetimer->timeIntervals[it]);CHKERRQ(ierr);        
	  }	  
	  thetimer->currInterval=0;
	  thetimer->numTimeSteps=thetimer->timeIntervals[thetimer->currInterval];
    }

    thetimer->startTimeStepResetFreq=-1;
    ierr = PetscOptionsGetInt(NULL,pre,"-start_time_step_reset_freq",&thetimer->startTimeStepResetFreq,&flg);CHKERRQ(ierr);
    if (flg) {
	  PetscInt tmp=0;
	  if (!thetimer->fixedStep) {
		for (it=0; it<thetimer->maxNumIntervals; it++) {
		  tmp=tmp+(thetimer->timeIntervals[it]);
        } 
	  } else {
        tmp=thetimer->numTimeSteps;
      }
      if (tmp > thetimer->startTimeStepResetFreq) {
        SETERRQ1(PETSC_COMM_WORLD,1,"Start time reset frequency less than total number of timer steps for StepTimer object %s",pre);
      } 
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time will be reset every %d steps for StepTimer object %s\n", thetimer->startTimeStepResetFreq, pre);CHKERRQ(ierr);	  
    }    
    
    thetimer->haveResetStartTimeStep=PETSC_FALSE;
	thetimer->count=0;
	    
    return 0;
    
}

PetscErrorCode updateStepTimer( const char pre[], PetscInt Iter, StepTimer *thetimer )
{

    PetscBool endOfSequence = PETSC_TRUE;
    PetscErrorCode ierr;
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Updating StepTimer object %s at iter %d\n", pre,Iter);CHKERRQ(ierr);        

	thetimer->count=0; /* reset counter */
    thetimer->haveResetStartTimeStep=PETSC_FALSE;
    
    if (!thetimer->fixedStep) {
      thetimer->currInterval++;
      if (thetimer->currInterval==thetimer->maxNumIntervals) {
//      We're now at the end of the sequence      
        thetimer->currInterval=0;
      } else {
//      Still within sequence      
        endOfSequence = PETSC_FALSE;
      }       
      thetimer->numTimeSteps=thetimer->timeIntervals[thetimer->currInterval];
      ierr = PetscPrintf(PETSC_COMM_WORLD,"New interval for StepTimer object %s at iter %d is %d\n", pre,Iter,thetimer->numTimeSteps);CHKERRQ(ierr);
    }

    if ((thetimer->startTimeStepResetFreq > 0) && (endOfSequence)) {
      thetimer->startTimeStep = thetimer->startTimeStep + thetimer->startTimeStepResetFreq;
      thetimer->haveResetStartTimeStep=PETSC_TRUE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"New start time step for StepTimer object %s is %d\n", pre, thetimer->startTimeStep);CHKERRQ(ierr);
    }

    return 0;
    
}

PetscErrorCode iniTimeDependentTimer( const char pre[], TimeDependentTimer *thetimer )
{

    PetscErrorCode ierr;
    PetscBool flg;
    PetscInt it;    
	PetscViewer fd;
	int fp;
    char timeFile[PETSC_MAX_PATH_LEN];

/*  read time data */
    ierr = PetscOptionsGetInt(NULL,pre,"-num_times",&thetimer->numTimes,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of time slices with the -%snum_times option",pre);

	ierr = PetscMalloc((thetimer->numTimes)*sizeof(PetscScalar), &thetimer->tdt);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,pre,"-times_file",timeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate name of time file with the -%stimes_file option",pre);
	
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,timeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&thetimer->tdt[0],thetimer->numTimes,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time-dependent object %s specified at times:\n",pre);CHKERRQ(ierr);            
	for (it=0; it<thetimer->numTimes; it++) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"tdt=%10.5f\n", thetimer->tdt[it]);CHKERRQ(ierr);        
	}    

    return 0;
    
}
