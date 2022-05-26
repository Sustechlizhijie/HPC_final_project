static char help[] = "Expilict Euler for 1D heat problem .\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <math.h>

#define pi acos(-1)   /* define pi */

int main(int argc,char **args)
{
  Vec            x, z, b;          
  Mat            A;                /* linear system matrix */
  PetscErrorCode ierr;
  PetscInt       i, n=200, start=0, end=n, col[3], rstart,rend,nlocal,rank;
  PetscReal      p=1.0, c=1.0, k=1.0, alpha, beta, dx, ix;
  PetscReal      dt=0.00001, t=0.0;
  PetscScalar    zero = 0.0, value[3] u0=0.0;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr);

  dx=1.0/n;
  alpha = k/p/c;
  beta = alpha*dt/dx/dx;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr); /* check the value */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"beta = %f\n",beta);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,n+1);CHKERRQ(ierr);
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
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1.0-2.0*beta; value[1] = beta;
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  if (rend == n+1) 
  {
    rend = n;
    i    = n; col[0] = n-1; col[1] = n; value[0] = beta; value[1] = 1.0-2.0*beta;
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
      for(int i=1; i<n; i++){
        ix = i*dx;
        u0 = exp(ix);
	      ierr = VecSetValues(z, 1, &i, &u0, INSERT_VALUES);CHKERRQ(ierr);
      }
  }
  
  ierr = VecAssemblyBegin(z);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(z);CHKERRQ(ierr);
  ix = 0.0;

  ierr = VecSet(b,zero);CHKERRQ(ierr);
  if(rank == 0){
    for(int i = 1; i < n+1; i++){
      PetscReal f;
      f = dt*sin(i*dx*pi);
      ierr = VecSetValues(b, 1, &i, &f, INSERT_VALUES);CHKERRQ(ierr);
    }
  }

  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  while(PetscAbsReal(t)<3.0){
     t += dt;
     ierr = MatMult(A,z,x);CHKERRQ(ierr);
     ierr = VecAXPY(x,1.0,b);CHKERRQ(ierr);

     ierr = VecSetValues(x, 1, &start, &zero, INSERT_VALUES);CHKERRQ(ierr);
     ierr = VecSetValues(x, 1, &end, &zero, INSERT_VALUES);CHKERRQ(ierr);
     ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
     ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
     ierr = VecCopy(x,z);CHKERRQ(ierr);

  }
  
  ierr = VecView(z,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
 
  ierr = VecDestroy(&x);CHKERRQ(ierr); 
  ierr = VecDestroy(&z);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();
  return ierr;
}

// EOF
