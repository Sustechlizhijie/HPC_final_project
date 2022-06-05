static char help[] = "impilict Euler for 1D heat problem .\n\n";

#include <petscksp.h>
#include <petscmath.h>
#include <petscsys.h>
#include <petscviewerhdf5.h>
#include <math.h>

#define pi acos(-1)   /* define pi */

int main(int argc,char **args)
{
  Vec            x, b, u, temp;          /* build the vecotr */
  Mat            A;                /* build the  matrix */
  KSP            ksp;
  PC             pc;
  PetscErrorCode ierr;             /* error checking */
  PetscInt       i, n=200, start=0, end=n, col[3], rstart,rend,nlocal,rank, iteration=0; /* n is region */
  PetscReal      p=1.0, c=1.0, k=1.0, alpha, beta, dx, ix, f;/* pck is the physic parameter */
  PetscReal      dt=0.00001, t=0.0, u0=0.0;   /* time step */
  PetscScalar    zero = 0.0, value[3], data[3];  /* u0 initial condition */
  PetscBool      restart = PETSC_FALSE;
  PetscViewer    h5;

  ierr = PetscInitialize(&argc,&args,(char*)0,help);if (ierr) return ierr;/* initial petsc */

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,NULL,NULL,NULL);CHKERRQ(ierr); 
  ierr = PetscOptionsGetReal(NULL,NULL,"-dt",&dt,NULL);CHKERRQ(ierr); /* read dt from command line */
  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr); /* read n from command line */
  ierr = PetscOptionsGetBool(NULL,NULL,"-restart",&restart,NULL);CHKERRQ(ierr); 
  ierr = PetscOptionsEnd();CHKERRQ(ierr); 

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);CHKERRQ(ierr); /* set up for MPI */
  ierr = PetscPrintf(PETSC_COMM_WORLD, "n = %d\n", n);CHKERRQ(ierr); /* print n */

  end=n; /* update the n value */

  dx=1.0/n;
  alpha = k/p/c;
  beta = alpha*dt/dx/dx;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dx = %f\n",dx);CHKERRQ(ierr); /* check the dx */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"dt = %f\n",dt);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"beta = %f\n",beta);CHKERRQ(ierr);/* check the beta */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"restart = %d\n",restart);CHKERRQ(ierr); 

  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);/* create vector */
  ierr = VecCreate(PETSC_COMM_WORLD,&temp);CHKERRQ(ierr); 

  ierr = VecSetSizes(x,PETSC_DECIDE,n+1);CHKERRQ(ierr); /* vector size*/
  ierr = VecSetSizes(temp, 3, PETSC_DECIDE);CHKERRQ(ierr);

  ierr = VecSetFromOptions(x);CHKERRQ(ierr);  /* enable option */
  ierr = VecSetFromOptions(temp);CHKERRQ(ierr);

  ierr = VecDuplicate(x,&b);CHKERRQ(ierr);  /* copy the type and layout */
  ierr = VecDuplicate(x,&u);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRQ(ierr); /* set start and end */
  ierr = VecGetLocalSize(x,&nlocal);CHKERRQ(ierr);  /* query the layout */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);  /* create matrix A */
  ierr = MatSetSizes(A,nlocal,nlocal,n+1,n+1);CHKERRQ(ierr); /* matrix size*/
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);  /* enable option */
  ierr = MatSetUp(A);CHKERRQ(ierr); 


  if (!rstart)      /* set the first line element */
  {
    rstart = 1;
    i      = 0; col[0] = 0; col[1] = 1; value[0] = 1.0+2.0*beta; value[1] = -beta; /* give values */
    ierr   = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr); /* set values */
  }
  
  if (rend == n+1)    /* set the final line element */
  {
    rend = n;
    i    = n; col[0] = n-1; col[1] = n; value[0] = -beta; value[1] = 1.0+2.0*beta;  /* give values */
    ierr = MatSetValues(A,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr); /* set values */
  }

  value[0] = -beta; value[1] = 1.0+2.0*beta; value[2] = -beta;   /* set the rest line element */
  for (i=rstart; i<rend; i++) 
  {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr   = MatSetValues(A,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr); /* set values */
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  /* Assemble the matrix */
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);   /* print the matrix */

  if(restart)
   {   
      ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"implicit_heat.h5", FILE_MODE_READ, &h5);CHKERRQ(ierr);    
      ierr = PetscObjectSetName((PetscObject) z, "implicit_heat_z");CHKERRQ(ierr);    
      ierr = PetscObjectSetName((PetscObject) temp, "implicit_heat_temp");CHKERRQ(ierr);
      ierr = VecLoad(temp, h5);CHKERRQ(ierr);    
      ierr = VecLoad(z, h5);CHKERRQ(ierr);    
      ierr = PetscViewerDestroy(&h5);CHKERRQ(ierr);  

      index=0;    
      ierr = VecGetValues(tem,1,&index,&dx);CHKERRQ(ierr);    
      index += 1;    
      ierr = VecGetValues(tem,1,&index,&dt);CHKERRQ(ierr);   
      index += 1;   
      ierr = VecGetValues(tem,1,&index,&t);CHKERRQ(ierr);   
      index= 0;   
    }
  else 
    {
      ierr = VecSet(b,zero);CHKERRQ(ierr);  /* set vector */
      if(rank == 0)
      {
          for(int i=1; i<n; i++){   /* from 1 to n-1 point*/
            u0 = exp(i*dx);  /* set u0 */
            ierr = VecSetValues(b, 1, &i, &u0, INSERT_VALUES);CHKERRQ(ierr);
          }
      }
      ierr = VecAssemblyBegin(b);CHKERRQ(ierr); /* Assemble the vector*/
      ierr = VecAssemblyEnd(b);CHKERRQ(ierr);
    }

  ierr = VecSet(u,zero);CHKERRQ(ierr); /* set initial vectot b */
  if(rank == 0){
    for(int i = 1; i < n; i++){  /* from 1 to n-1 point*/
      f = dt*sin(i*dx*pi); /* heat supply value */
      ierr = VecSetValues(u, 1, &i, &f, INSERT_VALUES);CHKERRQ(ierr); /* set value */
    }
  }

  ierr = VecAssemblyBegin(u);CHKERRQ(ierr);  /* Assemble the vector*/
  ierr = VecAssemblyEnd(u);CHKERRQ(ierr); 
  ierr = VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr); 


  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);    
  ierr = PCSetType(pc,PCJACOBI);CHKERRQ(ierr);   
  ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);    /*设置各种误差值*/
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);   



  while(PetscAbsReal(t)<3.0){   /* set the caculate time */
     t += dt;   /* time advance*/

     ierr = VecAXPY(b,1.0,u);CHKERRQ(ierr);    /*right hand size value*/
     ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);    /*sovle it*/

     ierr = VecSetValues(x, 1, &start, &zero, INSERT_VALUES);CHKERRQ(ierr); /* set value*/
     ierr = VecSetValues(x, 1, &end, &zero, INSERT_VALUES);CHKERRQ(ierr);
     ierr = VecAssemblyBegin(x);CHKERRQ(ierr); /* Assemble the vector*/
     ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
     ierr = VecCopy(x,b);CHKERRQ(ierr);  /* copy x into z*/

      iter += 1;    /*记录迭代次数*/
        if((iter%10)==0)
        {    /*如果迭代次数为10的倍数，即每迭代十次*/

          data[0] = dx; data[1] = dt; data[2] = t;    /*将值赋给数组*/
          ierr = VecSet(temp,zero);CHKERRQ(ierr);    /*初始化矩阵*/
          for(index=0;index<3;index++){    /*循环遍历数组，并将值赋给向量*/
            u0 = data[index];    /*将数组的值赋给u0*/
            ierr = VecSetValues(temp,1,&index,&u0,INSERT_VALUES);CHKERRQ(ierr);    /*将矩阵赋值给向量*/
          }
          ierr = VecAssemblyBegin(temp);CHKERRQ(ierr);    /*通知其余并行块将向量统一*/
          ierr = VecAssemblyEnd(temp);CHKERRQ(ierr);    /*结束通知*/

          ierr = PetscViewerCreate(PETSC_COMM_WORLD,&h5);CHKERRQ(ierr);    /*创建输出指针*/
          ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD,"implicit_heat.h5", FILE_MODE_WRITE, &h5);CHKERRQ(ierr);    /*创建输出文件*/
          ierr = PetscObjectSetName((PetscObject) b, "implicit_heat_z");CHKERRQ(ierr);    /*将z输出的名字命名为explicit-vector*/
          ierr = PetscObjectSetName((PetscObject) temp, "implicit_heat_temp");CHKERRQ(ierr);    /*将tem输出的名字命名为implicit-necess-data*/
          ierr = VecView(temp, h5);CHKERRQ(ierr);    /*tem输出到文件*/
          ierr = VecView(b, h5);CHKERRQ(ierr);    /*z输出到文件*/
          ierr = PetscViewerDestroy(&h5);CHKERRQ(ierr);    /*关闭输出*/
        }

  }

  ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);  /* view the z vector */
 
   /*Viewer to output in HDF5 format*/
  PetscViewer pv;
  PetscViewerCreate(PETSC_COMM_WORLD,&pv);
  PetscViewerASCIIOpen(PETSC_COMM_WORLD,"u_final_implicit.dat",&pv);
  VecView(b, pv);
  PetscViewerDestroy(&pv);
  /* deallocate the vector and matirx */
  ierr = VecDestroy(&temp);CHKERRQ(ierr); 
  ierr = VecDestroy(&x);CHKERRQ(ierr);  
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);

  ierr = PetscFinalize();  /* finish */
  return ierr;
}

// EOF
