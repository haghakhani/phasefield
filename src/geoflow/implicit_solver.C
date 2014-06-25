/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: implicit_solver.C 164 2013-06-18 15:27:22Z haghakhani $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "../header/hpfem.h"

struct ContData{

  HashTable*  El_Table;
  HashTable*  Node_Table;
  PetscInt*   Num_elem_proc;
  PetscInt    Total_elem,rank,size;
  TimeProps*  Timeptr;
  double      LapCoef,delta_t;

} myctx;

PetscErrorCode MakeRHS(ContData *myctx, Vec b);

#undef __FUNCT__
#define __FUNCT__ "implicit_solver" 

int implicit_solver(LaplacianData *Laplacian)
{
  Vec                x,b,xlocal;      /* approx solution, RHS */
  Mat                A;            /* linear system matrix */
  KSP                ksp;         /* KSP context */
  PetscReal          norm;//,val1,val2;         /* norm of solution error */
  PetscErrorCode     ierr;
  PetscInt           xsize,*num_elem_proc,*to,*from,its;
  PetscMPIInt        rank,size;
  KSPConvergedReason reason;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

  /* -------------------------------------------------------------------
     Compute the matrix and right-hand-side vector that define
     the linear system, Ax = b.
     ------------------------------------------------------------------- */

  ierr = PetscMalloc(size*sizeof(PetscInt),&num_elem_proc); CHKERRQ(ierr);
  num_elem_proc[rank]=num_nonzero_elem(Laplacian->El_Table);

  if (rank>0)
    MPI_Send(&num_elem_proc[rank], 1, MPI_INT, 0, 22, PETSC_COMM_WORLD);
  
  if (rank==0)
    for (int i=1;i<size;i++)
      MPI_Recv(&num_elem_proc[i], 1, MPI_INT, i, 22, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);

  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(num_elem_proc, size, MPI_INT, 0, PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);


  int total_elem=0,start_elem=0;

  for (int i=0;i<size;i++) 
    total_elem+=num_elem_proc[i];


  // we have to create a map between local and global vector
  myctx.El_Table     = Laplacian->El_Table;
  myctx.Node_Table   = Laplacian->NodeTable;
  myctx.Total_elem   = total_elem;
  myctx.Num_elem_proc= num_elem_proc;
  myctx.rank         = rank;
  myctx.size         = size;
  myctx.Timeptr      = Laplacian->timeprops;
  myctx.LapCoef      = Laplacian->LapCoef;
  myctx.delta_t      = Laplacian->delta_t;


  /*
     Create parallel vectors
   */
  ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
  ierr = VecSetType(b,VECSTANDARD); CHKERRQ(ierr);
  ierr = VecSetSizes(b,num_elem_proc[rank],total_elem);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);

  /*
     right-hand-side vector.
   */

  ierr = MakeRHS(&myctx, b);CHKERRQ(ierr);

  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = VecCopy(b,x);CHKERRQ(ierr);

  /* 
     Create and assemble parallel matrix
   */

  ierr = MatCreateShell(PETSC_COMM_WORLD,num_elem_proc[rank],num_elem_proc[rank],total_elem,total_elem, &myctx,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatLaplacian2D_Mult);CHKERRQ(ierr);

  /*
     Create linear solver context
   */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  /* 
     Set operators. Here the matrix that defines the linear system
     also serves as the preconditioning matrix.
   */
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPFGMRES);CHKERRQ(ierr);

  /*
     Set default preconditioner for this program to be block Jacobi.
     This choice can be overridden at runtime with the option
     -pc_type <type>
   */
  //  ierr = KSPSetTolerances(ksp,1.e-7,PETSC_DEFAULT,1.e9,3000);CHKERRQ(ierr);

  ierr = KSPSetTolerances(ksp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,5000);CHKERRQ(ierr);
  /* -------------------------------------------------------------------
     Solve the linear system
     ------------------------------------------------------------------- */
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  /*
     Solve the linear system
   */
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  /* -------------------------------------------------------------------
     Check solution and clean up
     ------------------------------------------------------------------- */

  /*
     Check the error
   */
  //ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
  if (rank==0){

  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp,&norm);CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Norm of error %G iterations %D\n",norm,its);CHKERRQ(ierr);
  //PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  ierr = KSPGetConvergedReason(ksp,&reason);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"kind of divergence is: ...........%D \n", reason);
  //PetscSynchronizedFlush(PETSC_COMM_WORLD);
  }
  /*

   */
  //ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  update_phi(Laplacian->El_Table,x);
  //ierr = VecRestoreArray(x, &xx);CHKERRQ(ierr);
  /* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
   */

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr); //ierr = PetscFree(phin);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr); //ierr = PCDestroy(&pc);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
  PetscFree(num_elem_proc);CHKERRQ(ierr);


  return 0;
}

//================================================================================
#undef __FUNCT__
#define __FUNCT__ "MatLaplacian2D_Mult"
/*
   Matrix-vector product subroutine for the 2D Laplacian.

 */
PetscErrorCode MatLaplacian2D_Mult(Mat A,Vec x,Vec y)
{
  PetscInt          elem,i,myid,*indices, *num_elem_proc,xsize;
  PetscScalar       *x_ptr,*y_ptr;
  PetscErrorCode    ierr;
  HashEntryPtr      currentPtr;
  Element           *Curr_El; 
  ContData           *myctx;  

  ierr = MatShellGetContext(A,(void**) &myctx);CHKERRQ(ierr);

  ierr = VecGetArray(x,&x_ptr);CHKERRQ(ierr);

  elem=0;
  HashEntryPtr      *buck = myctx->El_Table->getbucketptr();
  for(i=0; i< myctx->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  *(Curr_El->get_state_vars()) = x_ptr[elem];
	  elem++;
	}
	currentPtr=currentPtr->next;      	    
      }
    }


  ierr = VecRestoreArray(x,&x_ptr);CHKERRQ(ierr);

  MPI_Barrier(PETSC_COMM_WORLD);
  move_data(myctx->size, myctx->rank, myctx->El_Table, myctx->Node_Table, myctx->Timeptr);
  MPI_Barrier(PETSC_COMM_WORLD);

  //here we need comminucation between processors to transfer the information of the gohst elements and MPI_BARRIER

  elem=0;
  for(i=0; i<myctx->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  Curr_El->calc_lap_phi(myctx->El_Table,myctx->Node_Table);
	}
	currentPtr=currentPtr->next;      	    
      }
    }
  MPI_Barrier(PETSC_COMM_WORLD);

  ierr = VecGetArray(y,&y_ptr);CHKERRQ(ierr); 

  for(i=0; i<myctx->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  y_ptr[elem] = *(Curr_El->get_state_vars()) - (myctx->LapCoef) * (myctx->delta_t) * 
	    (*(Curr_El->get_lap_phi()) + *(Curr_El->get_lap_phi()+1));
	  elem++;
	}
	currentPtr=currentPtr->next;      	    
      }
    }
  ierr = VecRestoreArray(y,&y_ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
//=====================================================================================================================
int num_nonzero_elem(HashTable *El_Table){
  int               num=0,myid,i;
  HashEntryPtr      currentPtr;
  Element           *Curr_El; 
  HashEntryPtr      *buck = El_Table->getbucketptr();
  //printf("the numebr of buckets are: %d", El_Table->get_no_of_buckets());

  for(i=0;i<El_Table->get_no_of_buckets();i++)
    if(*(buck+i)){
      currentPtr = *(buck+i);
      while(currentPtr){
	Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  num++;	
	}
	currentPtr=currentPtr->next;      	    
      }
    }

  return (num);
}
//==========================================================================================
#undef __FUNCT__
#define __FUNCT__ "MakeRHS"
/*
   Making the RHS of the Ax=b
 */
PetscErrorCode MakeRHS(ContData *ctx,Vec b){

  //PetscFunctionBegin;

  int            elem=0;
  HashTable      *El_Table=ctx->El_Table;
  HashEntryPtr   *buck = ctx->El_Table->getbucketptr();
  HashEntryPtr   currentPtr;
  Element        *Curr_El; 
  PetscInt       *num_elem_proc,rank,size,xsize;
  PetscScalar    *b_ptr;
  VecScatter     *vscat;
  Vec            blocal;
  PetscErrorCode ierr;

  ierr = VecGetArray(b,&b_ptr);CHKERRQ(ierr);

  for(int i=0; i<ctx->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  b_ptr[elem]=*(Curr_El->get_state_vars());
	  elem++;
	}

	currentPtr=currentPtr->next;      	    
      }
    }
  ierr = VecRestoreArray(b,&b_ptr);CHKERRQ(ierr);

  PetscFunctionReturn(0);

}

//==========================================================================================
PetscErrorCode update_phi(HashTable *El_Table, Vec update){

  int elem=0,myid,i;
  HashEntryPtr   *buck = El_Table->getbucketptr();
  HashEntryPtr   currentPtr;
  Element        *Curr_El; 
  PetscScalar    *update_ptr;
  PetscErrorCode ierr;


  ierr = VecGetArray(update,&update_ptr);CHKERRQ(ierr);

  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  *(Curr_El->get_state_vars())=update_ptr[elem];
	  elem++;
	}
	currentPtr=currentPtr->next;      	    
      }
    }
  ierr = VecRestoreArray(update,&update_ptr);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}
