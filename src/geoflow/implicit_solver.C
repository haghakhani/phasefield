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


#undef __FUNCT__
#define __FUNCT__ "implicit_solver" 

int implicit_solver(LaplacianData *Laplacian)
{
  Vec            x,b,w;      /* approx solution, RHS */
  Mat            A;            /* linear system matrix */
  KSP            ksp;         /* KSP context */
  //PC             pc;           /* PC context */
  PetscReal      norm,val1,val2;         /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       num_elem,its;
  PetscMPIInt    rank,size;
  PetscScalar    *phin,*xx,sizo;
  KSPConvergedReason reason;
     
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

  /* -------------------------------------------------------------------
     Compute the matrix and right-hand-side vector that define
     the linear system, Ax = b.
     ------------------------------------------------------------------- */

  num_elem=num_nonzero_elem(Laplacian->El_Table);
  //printf("Number of elements are (hi i am second)...........%d\n", num_nonzero_elem(Laplacian->El_Table));
   
  /*
    Create parallel vectors
  */
  ierr = VecCreate(MPI_COMM_SELF,&b);CHKERRQ(ierr);
  ierr = VecSetType(b,VECSEQ); CHKERRQ(ierr);
  ierr = VecSetSizes(b,PETSC_DECIDE,num_elem);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);
  /*
    right-hand-side vector.
  */
  
  ierr = MakeRHS(Laplacian,num_elem, b);CHKERRQ(ierr);
  //phin = phase(Laplacian->El_Table,num_elem);
  //ierr = VecPlaceArray(b,phin);CHKERRQ(ierr);

  ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
  ierr = VecCopy(b,x);CHKERRQ(ierr);

  // ierr = VecDuplicate(b,&w);CHKERRQ(ierr);
  // ierr = VecWAXPY(w,-1,x,b);CHKERRQ(ierr);

  // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"It is w \n");CHKERRQ(ierr);
  // ierr = VecView(w,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

  //  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"It is x \n");CHKERRQ(ierr);
  //  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  //  ierr = VecView(x,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);


  //  ierr = VecCreate(MPI_COMM_SELF,&x);CHKERRQ(ierr);
  //  ierr = VecSetType(x,VECSEQ); CHKERRQ(ierr);
  //  ierr = VecSetSizes(x,PETSC_DECIDE,num_elem);CHKERRQ(ierr);
  //  ierr = VecSet(x,1.0);

  //  PetscInt sizi[4];
  /* 
     Create and assemble parallel matrix
  */

  ierr = MatCreateShell(MPI_COMM_SELF,num_elem,num_elem,PETSC_DETERMINE,PETSC_DETERMINE, Laplacian,&A);CHKERRQ(ierr);
  ierr = MatShellSetOperation(A,MATOP_MULT,(void(*)(void))MatLaplacian2D_Mult);CHKERRQ(ierr);

  /*
    Create linear solver context
  */
  ierr = KSPCreate(MPI_COMM_SELF,&ksp);CHKERRQ(ierr);

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
  
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp,&norm);CHKERRQ(ierr);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Norm of error %G iterations %D\n",norm,its);CHKERRQ(ierr);
  PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  ierr = KSPGetConvergedReason(ksp,&reason);
  ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"kind of divergence is: ...........%D \n", reason);
  PetscSynchronizedFlush(PETSC_COMM_WORLD);
  /*
  
   */
  ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  update_phi(Laplacian->El_Table,xx);
  ierr = VecRestoreArray(x, &xx);CHKERRQ(ierr);
  /* 
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */

  ierr = KSPDestroy(&ksp);CHKERRQ(ierr); //ierr = PetscFree(phin);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr); //ierr = PCDestroy(&pc);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr); ierr = MatDestroy(&A);CHKERRQ(ierr);
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
  LaplacianData     *Laplacian;
  PetscInt          elem,i,myid,*indices, num_elem;
  PetscScalar       *xx,*yy;
  PetscErrorCode    ierr;
  HashEntryPtr      currentPtr;
  Element           *Curr_El; 
 
  PetscFunctionBegin;

  //   ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"Hello \n");CHKERRQ(ierr);
  //   PetscSynchronizedFlush(PETSC_COMM_WORLD); 

  ierr = MatShellGetContext(A,/*(void**)&*/&Laplacian);CHKERRQ(ierr);

  // ierr =  VecGetSize( x,&i);CHKERRQ(ierr);
  // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"this is x size= %d \n",i);CHKERRQ(ierr);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD); 

  // ierr =  VecGetSize( y,&i);CHKERRQ(ierr);
  // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"this is y size = %d\n",i);CHKERRQ(ierr);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD); 

  //ierr = VecView(x,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

  num_elem=num_nonzero_elem(Laplacian->El_Table);
  // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"this is number of element = %d\n",num_elem);CHKERRQ(ierr);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD); 

  PetscMalloc(num_elem*sizeof(PetscInt),&indices);
  PetscMalloc(num_elem*sizeof(PetscScalar),&xx);
  PetscMalloc(num_elem*sizeof(PetscScalar),&yy);

  for(i=0;i<num_elem;i++)
    indices[i]=i;

  //  ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  //  ierr = VecGetArray(y,&yy);CHKERRQ(ierr);
  //  ierr = VecAssemblyBegin(x);
  //  ierr = VecAssemblyEnd(x);
  ierr = VecGetValues(x,num_elem,indices,xx);CHKERRQ(ierr);
  elem=0;
  HashEntryPtr      *buck = Laplacian->El_Table->getbucketptr();
  for(i=0; i<Laplacian->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
     while(currentPtr){
  	Element* Curr_El=(Element*)(currentPtr->value);
  	if(Curr_El->get_adapted_flag()>0){
	  *(Curr_El->get_state_vars()) = xx[elem];
//  	  for(int j=0;j<4;j++)
//  	    if(*(Curr_El->get_neigh_proc()+j) == INIT)   // this is a boundary!
//  	      *(Curr_El->get_state_vars()) = -1;
          //if (*(Curr_El->get_state_vars()) != xx[elem]) 
	  //printf("this is in mult not equal xx[%d]=%f   phi=%f\n",  elem, xx[elem], *(Curr_El->get_state_vars()));
	  elem++;
  	}
  	currentPtr=currentPtr->next;      	    
      }
    }

  // ierr = VecGetArray(x,&xx);CHKERRQ(ierr);
  // ierr = VecGetArray(y,&yy);CHKERRQ(ierr);
  //int n = sizeof(a) / sizeof(int);
  //  elem=0;
  // for(i=0; i<Laplacian->El_Table->get_no_of_buckets();i++)
  //   if(*(buck+i)){
  //     currentPtr = *(buck+i);
  //     while(currentPtr){
  // 	Element* Curr_El=(Element*)(currentPtr->value);
  // 	if(Curr_El->get_adapted_flag()>0){
  // 	  //if (*(Curr_El->get_state_vars()) != xx[elem]) printf("not equal xx[elem]=%f   phi=%f\n",  xx[elem], *(Curr_El->get_state_vars()));
  // 	  *(Curr_El->get_state_vars()) = xx[elem];
  // 	  //if(xx[elem]>1) printf("x[elem]=%f \n",xx[elem]);
  // 	  for(int j=0;j<4;j++)
  // 	    if(*(Curr_El->get_neigh_proc()+j) == INIT)   // this is a boundary!
  // 	      *(Curr_El->get_state_vars()) = 0;
  // 	  elem++;
  // 	}    
  // 	currentPtr=currentPtr->next;      	    
  //     }
  //   }
  
  elem=0;
  for(i=0; i<Laplacian->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  Curr_El->calc_lap_phi(Laplacian->El_Table,Laplacian->NodeTable);
	}
	currentPtr=currentPtr->next;      	    
      }
    }
  for(i=0; i<Laplacian->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
  	Element* Curr_El=(Element*)(currentPtr->value);
  	if(Curr_El->get_adapted_flag()>0){
  	  yy[elem] = *(Curr_El->get_state_vars()) - (Laplacian->LapCoef) * (Laplacian->delta_t) * 
  	    (*(Curr_El->get_lap_phi()) + *(Curr_El->get_lap_phi()+1));
	  //	    if ( *(Curr_El->get_state_vars())>0)
	  //	    printf("phi=%f, lap_phi_x=%f, lap_phi_y=%f,  y=%f \n ", 
	  //	   *(Curr_El->get_state_vars()), *(Curr_El->get_lap_phi()),*(Curr_El->get_lap_phi()+1), yy[elem]);
	  elem++;
  	}
  	currentPtr=currentPtr->next;      	    
      }
    }
  ierr=VecSetValues(y,num_elem,indices,yy,INSERT_VALUES);
  VecAssemblyBegin(y);
  VecAssemblyEnd(y);
  // ierr = VecRestoreArray(x, &xx);CHKERRQ(ierr);
  // ierr = VecRestoreArray(y, &yy);CHKERRQ(ierr);

  PetscFree(xx);
  PetscFree(yy);
  PetscFree(indices);
  //exit(1); 
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
PetscErrorCode MakeRHS(LaplacianData *Laplacian,PetscInt num_elem,Vec b){
  PetscFunctionBegin;

  int elem=0,myid,i;
  // HashTable *El_Table;
  HashEntryPtr      *buck = Laplacian->El_Table->getbucketptr();
  HashEntryPtr      currentPtr;
  Element           *Curr_El; 
  PetscScalar       *phin,*xx;
  PetscInt       *indices;
  PetscErrorCode  ierr;
 

  PetscMalloc(num_elem*sizeof(PetscScalar),&phin);
  PetscMalloc(num_elem*sizeof(PetscScalar),&xx);
  PetscMalloc(num_elem*sizeof(PetscInt),&indices);
  for(i=0;i<num_elem;i++)
    indices[i]=i;
  
  for(i=0; i<Laplacian->El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  phin[elem]=*(Curr_El->get_state_vars());
	  //if ( phin[elem]!=*(Curr_El->get_state_vars())) printf("phi[%d]=%f    state_vars=%f \n",elem,phin[elem],*(Curr_El->get_state_vars()));
	  elem++;
	}
   
	currentPtr=currentPtr->next;      	    
      }
    }

  ierr=VecSetValues(b,num_elem,indices,phin,INSERT_VALUES);
  VecAssemblyBegin(b);
  VecAssemblyEnd(b);
  //=======================================
  // ierr = VecGetValues(b,num_elem,indices,xx);CHKERRQ(ierr);
  // elem=0;
  // //  HashEntryPtr      *buck = Laplacian->El_Table->getbucketptr();
  // for(i=0; i<Laplacian->El_Table->get_no_of_buckets(); i++)
  //   if(*(buck+i)){
  //     currentPtr = *(buck+i);
  //     while(currentPtr){
  // 	Curr_El=(Element*)(currentPtr->value);
  // 	if(Curr_El->get_adapted_flag()>0){
  //         if (/**(Curr_El->get_state_vars()) != xx[elem] && */elem==266) 
  // 	    printf("this is in RHS not equal xx[%d]=%f   phi=%f\n",  elem, xx[elem], *(Curr_El->get_state_vars()));
  // 	  elem++;
  // 	}
  // 	currentPtr=currentPtr->next;      	    
  //     }
  //   }
  // ierr = PetscSynchronizedPrintf(PETSC_COMM_WORLD,"It is b \n");CHKERRQ(ierr);
  // PetscSynchronizedFlush(PETSC_COMM_WORLD); 
  // ierr = VecView(b,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

  PetscFree(xx);
  PetscFree(phin);
  PetscFree(indices);
  PetscFunctionReturn(0);

}

//=========================================================================================
// PetscScalar *phase(HashTable *El_Table, int num_elem){
//   int elem=0,myid,i;
//   HashEntryPtr      *buck = El_Table->getbucketptr();
//   HashEntryPtr      currentPtr;
//   Element           *Curr_El; 
//   PetscScalar       *phin;

//   PetscMalloc(num_elem*sizeof(PetscScalar),&phin);
  
//   for(i=0; i<El_Table->get_no_of_buckets(); i++)
//     if(*(buck+i)){
//       HashEntryPtr currentPtr = *(buck+i);
//       while(currentPtr){
// 	Element* Curr_El=(Element*)(currentPtr->value);
// 	if(Curr_El->get_adapted_flag()>0){
// 	  phin[elem]=*(Curr_El->get_state_vars());
// 	  elem++;
// 	}
   
// 	currentPtr=currentPtr->next;      	    
//       }
//     }

//   PetscFunctionReturn(phin);
// }

//==========================================================================================
void update_phi(HashTable *El_Table, double *update){
  int elem=0,myid,i;
  HashEntryPtr      *buck = El_Table->getbucketptr();
  HashEntryPtr      currentPtr;
  Element           *Curr_El; 

  for(i=0; i<El_Table->get_no_of_buckets(); i++)
    if(*(buck+i)){
      HashEntryPtr currentPtr = *(buck+i);
      while(currentPtr){
	Element* Curr_El=(Element*)(currentPtr->value);
	if(Curr_El->get_adapted_flag()>0){
	  *(Curr_El->get_state_vars())=update[elem];
	  elem++;
	}
	currentPtr=currentPtr->next;      	    
      }
    }

  return;
}