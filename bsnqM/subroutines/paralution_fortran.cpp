// **************************************************************************
//
//    PARALUTION   www.paralution.com
//
//    Copyright (C) 2015  PARALUTION Labs UG (haftungsbeschränkt) & Co. KG
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRA 706051
//                        Vertreten durch:
//                        PARALUTION Labs Verwaltungs UG (haftungsbeschränkt)
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRB 721277
//                        Geschäftsführer: Dimitar Lukarski, Nico Trost
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// **************************************************************************



// PARALUTION version 1.1.0 


// #######################################################################################################
// ###                                                                                                 ###
// ###                               PARALUTION FORTRAN PLUG-IN                                        ###
// ###                                                                                                 ###
// ###                                                                                                 ###
// ###     Each function listed here can be executed from Fortran code by simply calling it:           ###
// ###                                                                                                 ###
// ###        C function name: "void paralution_fortran_function_();"                                  ###
// ###        Fortran syntax: "call paralution_fortran_function()"                                     ###
// ###                                                                                                 ###
// #######################################################################################################
// ###                                                                                                 ###
// ###     paralution_init                                                                             ###
// ###        Initalize the PARALUTION backend                                                         ###
// ###                                                                                                 ###
// #######################################################################################################
// ###                                                                                                 ###
// ###     paralution_stop                                                                             ###
// ###        Stops the PARALUTION backend                                                             ###
// ###                                                                                                 ###
// #######################################################################################################
// ###                                                                                                 ###
// ###     paralution_fortran_solve_coo_                                                               ###
// ###        Solves a linear system for given COO matrix, rhs, solution vector,                       ###
// ###        solver and preconditioner.                                                               ###
// ###                                                                                                 ###
// ###        input parameters:                                                                        ###
// ###                                                                                                 ###
// ###          n          -  number of matrix rows                                                    ###
// ###          m          -  number of matrix cols                                                    ###
// ###          nnz        -  number of non-zero matrix elements                                       ###
// ###          solver     -  solver string, can be CG,BiCGStab,FixedPoint,GMRES,FGMRES                ###
// ###          mformat    -  matrix format for the solving procedure, can be                          ###
// ###                        CSR,MCSR,BCSR,COO,DIA,ELL,HYB,DENSE                                      ###
// ###          precond    -  preconditioner string, can be None, Jacobi, MultiColoredGS,              ###
// ###                        MultiColoredSGS, ILU, MultiColoredILU, FSAI                              ###
// ###          pformat    -  preconditioner format for MultiColored preconditioners, can be           ###
// ###                        CSR,MCSR,BCSR,COO,DIA,ELL,HYB,DENSE                                      ###
// ###          row        -  matrix row index - see COO format                                        ###
// ###          col        -  matrix col index - see COO format                                        ###
// ###          val        -  matrix values    - see COO format                                        ###
// ###          rhs        -  right hand side vector                                                   ###
// ###          atol       -  absolute tolerance                                                       ###
// ###          rtol       -  relative tolerance                                                       ###
// ###          maxiter    -  maximum iterations allowed                                               ###
// ###          basis      -  Basis size when using GMRES                                              ###
// ###          p          -  ILU(p) factorization based on power                                      ###
// ###                        p needs to be specified when ILU preconditioner is chosen.               ###
// ###          q          -  ILU(p,q)                                                                 ###
// ###                        p and q need to be specified when MultiColoredILU preconditioner         ###
// ###                        is chosen.                                                               ###
// ###                                                                                                 ###
// ###        output parameters:                                                                       ###
// ###                                                                                                 ###
// ###          x          -  solution vector                                                          ###
// ###          iter       -  iteration count                                                          ###
// ###          resnorm    -  residual norm                                                            ###
// ###          err        -  error code                                                               ###
// ###                                                                                                 ###
// #######################################################################################################
// ###                                                                                                 ###
// ###     paralution_fortran_solve_csr_                                                               ###
// ###        Solves a linear system for given CSR matrix, rhs, solution vector,                       ###
// ###        solver and preconditioner.                                                               ###
// ###                                                                                                 ###
// ###        input parameters:                                                                        ###
// ###                                                                                                 ###
// ###          n          -  number of matrix rows                                                    ###
// ###          m          -  number of matrix cols                                                    ###
// ###          nnz        -  number of non-zero matrix elements                                       ###
// ###          solver     -  solver string, can be CG,BiCGStab,FixedPoint,GMRES                       ###
// ###          mformat    -  matrix format for the solving procedure, can be                          ###
// ###                        CSR,MCSR,BCSR,COO,DIA,ELL,HYB,DENSE                                      ###
// ###          precond    -  preconditioner string, can be None, Jacobi, MultiColoredGS,              ###
// ###                        MultiColoredSGS, ILU, MultiColoredILU                                    ###
// ###          pformat    -  preconditioner format for MultiColored preconditioners, can be           ###
// ###                        CSR,MCSR,BCSR,COO,DIA,ELL,HYB,DENSE                                      ###
// ###          row_offset -  matrix row offset - see CSR format                                       ###
// ###          col        -  matrix col index  - see CSR format                                       ###
// ###          val        -  matrix values     - see CSR format                                       ###
// ###          rhs        -  right hand side vector                                                   ###
// ###          atol       -  absolute tolerance                                                       ###
// ###          rtol       -  relative tolerance                                                       ###
// ###          maxiter    -  maximum iterations allowed                                               ###
// ###          basis      -  Basis size when using GMRES                                              ###
// ###          p          -  ILU(p) factorization based on power                                      ###
// ###                        p needs to be specified when ILU preconditioner is chosen.               ###
// ###          q          -  ILU(p,q)                                                                 ###
// ###                        p and q need to be specified when MultiColoredILU preconditioner         ###
// ###                        is chosen.                                                               ###
// ###                                                                                                 ###
// ###        output parameters:                                                                       ###
// ###                                                                                                 ###
// ###          x          -  solution vector                                                          ###
// ###          iter       -  iteration count                                                          ###
// ###          resnorm    -  residual norm                                                            ###
// ###          err        -  error code                                                               ###
// ###                                                                                                 ###
// #######################################################################################################
// ###                                                                                                 ###
// ###        error codes:                                                                             ###
// ###                                                                                                 ###
// ###           0  -  no error                                                                        ###
// ###           1  -  absolute tolerance is reached                                                   ###
// ###           2  -  relative tolerance is reached                                                   ###
// ###           3  -  divergence tolerance is reached                                                 ###
// ###           4  -  max iter is reached                                                             ###
// ###           5  -  invalid solver                                                                  ###
// ###           6  -  invalid preconditioner                                                          ###
// ###           7  -  invalid matrix format                                                           ###
// ###           8  -  invalid preconditioner format                                                   ###
// ###                                                                                                 ###
// #######################################################################################################
#include <paralution.hpp>


//-------------------------------Types-------------------------------
class lsClass
{
public:

  int     np, nnz, t=0;
  int     maxiter;
  //CSR Storage
  // int     *iv, *jv;
  // double  *mat;  
  double  *rhs, *solX;
  double  atol, rtol, div;

  paralution::IterativeLinearSolver<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > *ls = NULL;

  // paralution::LocalVector<double> par_x;
  // paralution::LocalVector<double> par_rhs;
  paralution::LocalMatrix<double> par_mat;
  

  lsClass()
  {
    np = 0;
    nnz = 0;
    // iv = NULL;
    // jv = NULL;
    // mat = NULL;
    rhs = NULL;
    solX = NULL;
    return;
  }

  void init(int, int, const int*, const int*, const double*,
    double, double, double, int );
  void check();
  void solve(const double*, double*, int&, double&, int&);  

};



typedef void *opqObj;

//-----------------------------End Types-----------------------------



//------------------------------Declare------------------------------
extern "C" 
{
  // Shagun modify 2017_08_14
  void paralution_init(int);
  void paralution_stop(void);  
  void paralution_fortran_solve_csr(int, int, int, char*, char*, char*, char*, const int*, const int*,
                                    const double*, const double*, double, double, double, int, int,
                                    int, int, double*, int&, double&, int&);

  opqObj  createLSObj();
  void    initLSObj(opqObj, int, int, const int*, const int*, 
    const double*, double, double, double, int );
  void    checkLSObj(opqObj);
  void    solveLSObj(opqObj, const double*, double*, int&, double&, int&);
}

void paralution_fortran_solve(double, double, double, int,
                              paralution::LocalMatrix<double>*, paralution::LocalVector<double>*,
                              paralution::LocalVector<double>*, int&, double&, int&);

//----------------------------End Declare----------------------------



//------------------------------opqObj-------------------------------
opqObj createLSObj()
{
  lsClass *lsObj = new lsClass();
  return (opqObj)lsObj;
}



void initLSObj(opqObj vObj, int np, int nnz, const int *f90_iv, 
  const int *f90_jv, const double *f90_mat, double in_atol, 
  double in_rtol, double in_div, int in_maxiter)
{
  lsClass *lsObj = (lsClass *)vObj;

  lsObj->init(np, nnz, f90_iv, f90_jv, f90_mat, 
    in_atol, in_rtol, in_div, in_maxiter);
}



void checkLSObj(opqObj vObj)
{
  lsClass *lsObj = (lsClass *)vObj;

  lsObj->check();
}



void solveLSObj(opqObj vObj, const double *f90_rhs, double *f90_x, 
  int &iter, double &resnorm, int &err)
{

  lsClass *lsObj = (lsClass *)vObj;

  lsObj->solve(f90_rhs, f90_x,
    iter, resnorm, err);

}            
//----------------------------End opqObj-----------------------------



//------------------------------lsClass------------------------------
void lsClass::init(int in_np, int in_nnz, const int *in_iv, 
  const int *in_jv, const double *in_mat, double in_atol, 
  double in_rtol, double in_div, int in_maxiter)
{

  // paralution::LocalMatrix<double> par_mat;  

  int *iv     = NULL;
  int *jv     = NULL;
  double *mat = NULL;

  np = in_np;
  nnz = in_nnz;

  atol = in_atol;
  rtol = in_rtol;
  div  = in_div;
  maxiter = in_maxiter;

  paralution::allocate_host(np+1, &iv);
  paralution::allocate_host(nnz, &jv);
  paralution::allocate_host(nnz, &mat);

  paralution::allocate_host(np, &rhs);
  paralution::allocate_host(np, &solX);

  // Shift since Fortran arrays start at 1
  for (int i=0; i<np+1; ++i)  
    iv[i] = in_iv[i] - 1;

  for (int i=0; i<nnz; ++i) {
    jv[i] = in_jv[i] - 1;
    mat[i] = in_mat[i];
  }  

  // Allocate paralution data structures
  par_mat.SetDataPtrCSR(&iv, &jv, &mat, "Imported Fortran CSR Matrix", nnz, np, np);
  par_mat.info();

  ls = new paralution::BiCGStab<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double >;

  ls->SetOperator(par_mat);
  ls->Init(atol, rtol, div, maxiter);    
  ls->Build();

  // delete [] iv, jv, mat;
  return;
}



void lsClass::check()
{
  t++;
  // std::cout << iv[np] << " " << t <<"\n";
  std::cout << np << " " << t <<"\n";
  return;  
}



void lsClass::solve(const double *in_rhs, double *in_x, 
  int &iter, double &resnorm, int &err)
{
  paralution::LocalVector<double> par_x;
  paralution::LocalVector<double> par_rhs;
  // paralution::LocalMatrix<double> par_mat; 

  for (int i=0; i<np; ++i)
    rhs[i] = in_rhs[i];

  for (int i=0; i<np; ++i)  
    solX[i] = in_x[i];

  par_rhs.SetDataPtr(&rhs, "Imported Fortran rhs", np);    
  par_x.SetDataPtr(&solX, "Imported Fortran x", np);

  // // Allocate paralution data structures
  // par_mat.SetDataPtrCSR(&iv, &jv, &mat, "Imported Fortran CSR Matrix", nnz, np, np);
  // par_mat.info();

  // paralution_fortran_solve(atol, rtol, div, maxiter,
  //                          &par_mat, &par_rhs, &par_x, iter, resnorm, err);

  // ls->SetOperator(par_mat);
  // ls->Init(atol, rtol, div, maxiter);  
  // ls->Build();

  // par_mat->ConvertToCSR();  

  ls->Solve(par_rhs, &par_x);

  iter = ls->GetIterationCount();
  resnorm = ls->GetCurrentResidual();
  err = ls->GetSolverStatus();

  // ls->Clear();

  par_x.MoveToHost();
  par_x.LeaveDataPtr(&solX);

  par_rhs.LeaveDataPtr(&rhs);
  // par_mat.LeaveDataPtrCSR(&iv, &jv, &mat);


  for (int i=0; i<np; ++i)
    in_x[i] = solX[i];

  return;
}
//----------------------------End lsClass----------------------------



//----------------------------Paralution-----------------------------
// Shagun modify 2017_08_14
/// Initializes the PARALUTION backend
void paralution_init(int nthreads) {

  paralution::init_paralution();
  paralution::set_omp_threads_paralution(nthreads);  
  paralution::info_paralution();

}



/// Stops the PARALUTION backend
void paralution_stop(void) {

  paralution::stop_paralution();

}



/// Solves a linear system for given CSR matrix, rhs, solution vector, solver and preconditioner.
void paralution_fortran_solve_csr(int n, int m, int nnz, char *solver, char *mformat, char *precond, char *pformat,
                                  const int *fortran_row_offset, const int *fortran_col, const double *fortran_val,
                                  const double *fortran_rhs, double atol, double rtol, double div, int maxiter,
                                  int basis, int p, int q, double *fortran_x, int &iter, double &resnorm, int &err) {

  paralution::LocalVector<double> paralution_x;
  paralution::LocalVector<double> paralution_rhs;
  paralution::LocalMatrix<double> paralution_mat;

  int *row_offset = NULL;
  int *col        = NULL;
  double *val     = NULL;

  paralution::allocate_host(n+1, &row_offset);
  paralution::allocate_host(nnz, &col);
  paralution::allocate_host(nnz, &val);

  double *in_rhs = NULL;
  double *in_x = NULL;
  paralution::allocate_host(m, &in_rhs);
  paralution::allocate_host(n, &in_x);

  for (int i=0; i<m; ++i)
    in_rhs[i] = fortran_rhs[i];
  for (int i=0; i<n; ++i)
    in_x[i] = fortran_x[i];

  paralution_rhs.SetDataPtr(&in_rhs, "Imported Fortran rhs", m);
  paralution_x.SetDataPtr(&in_x, "Imported Fortran x", n);

  // Copy matrix so we can convert it to any other format without breaking the fortran code
  // Shift since Fortran arrays start at 1
  for (int i=0; i<n+1; ++i)
    row_offset[i] = fortran_row_offset[i] - 1;

  for (int i=0; i<nnz; ++i) {
    col[i] = fortran_col[i] - 1;
    val[i] = fortran_val[i];
  }

  // Allocate paralution data structures
  paralution_mat.SetDataPtrCSR(&row_offset, &col, &val, "Imported Fortran CSR Matrix", nnz, n, m);
  paralution_mat.info();

  paralution_fortran_solve(atol, rtol, div, maxiter,
                           &paralution_mat, &paralution_rhs, &paralution_x, iter, resnorm, err);

  paralution_x.MoveToHost();
  paralution_x.LeaveDataPtr(&in_x);

  for (int i=0; i<n; ++i)
    fortran_x[i] = in_x[i];

  delete [] in_x;

}


void paralution_fortran_solve(double atol, double rtol,
                              double div, int maxiter, paralution::LocalMatrix<double> *mat,
                              paralution::LocalVector<double> *rhs, paralution::LocalVector<double> *x,
                              int &iter, double &resnorm, int &err) {

  // Iterative Linear Solver and Preconditioner
  paralution::IterativeLinearSolver<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double > *ls = NULL;  
  
  // Switch for solver selection
  ls = new paralution::BiCGStab<paralution::LocalMatrix<double>, paralution::LocalVector<double>, double >;
  

  // Switch for preconditioner selection
  

  ls->SetOperator(*mat);
  ls->Init(atol, rtol, div, maxiter);  

  ls->Build();

  mat->ConvertToCSR();  

  x->info();
  rhs->info();
  mat->info();

  ls->Solve(*rhs, x);

  iter = ls->GetIterationCount();
  resnorm = ls->GetCurrentResidual();
  err = ls->GetSolverStatus();

  ls->Clear();
  delete ls;  
}
//--------------------------End Paralution---------------------------