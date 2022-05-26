static char help[] = "Expilict Euler for 1D heat problem .\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <math.h>



int main(int argc,char **args)
{
  Vec            x, z, b;          
  Mat            A;                /* linear system matrix */
  PetscReal      norm = 0.0, normt = 1.0, tol=1000.*PETSC_MACHINE_EPSILON;  /* norm of solution error */
  PetscErrorCode ierr;
  PetscInt       i,n = 201 ,col[3],rstart,rend,nlocal,rank;
  PetscReal      dx=1/(n-1), dt=0.00001, p=1.0, c=1.0, k=1.0;
  PetscScalar    zero = 0.0, value[3];

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

  alpha = k/p/c
  beta = alpha*dt/dx/dx

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&z);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr);
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,nlocal,nlocal,n,n);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);


  if (!rstart) 
  {
    rstart = 1;
    i      = 0; col[0] = 0.0; col[1] = 1.0; value[0] = 1.0-2.0*beta; value[1] = beta;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n) 
  {
    rend = n-1;
    i    = n-1; col[0] = n-2; col[1] = n-1; value[0] = beta; value[1] = 1.0-2.0*beta;
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  value[0] = beta; value[1] = 1.0-2.0*beta; value[2] = beta;
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }

  /* Assemble the matrix */
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


  ierr = VecSet(z,zero);CHKERRQ(ierr);
  if(rank == 0){
      for(int i=0; i<n; i++){
        PetscReal u0;
        i = 0;
	    ierr = VecSetValues(z, 1, &i, &one, INSERT_VALUES);CHKERRQ(ierr);
      }
	  
  }
  
  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);
  
  ierr = VecSet(b,zero);CHKERRQ(ierr);
  if(rank == 0){
    for(int i = 0; i < n; i++){
	  PetscReal inp;
      inp = dt*sin(pi*(i+1.0)*dx);
	  ierr = VecSetValues(b, 1, &i, &inp, INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  while(PetscAbsReal(norm-normt)>tol){
     normt= norm;
     ierr = MatMult(A,z,x);CHKERRQ(ierr);
     ierr = VecNorm(x,NORM_2,&norm);CHKERRQ(ierr);
     ierr = VecScale(x,(PetscScalar)1.0/norm);CHKERRQ(ierr);
     ierr = VecAXPY(x,1.0,b);CHKERRQ(ierr);
     ierr = VecCopy(x,z);CHKERRQ(ierr);
  }
  
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"solution = %f\n",norm);CHKERRQ(ierr);
 
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

// EOF
