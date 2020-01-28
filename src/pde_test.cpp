/* =========================================================================
   Copyright (c) 2016-2020, Department of Engineering Physics,
                            Tsinghua University, Beijing, China.

   Portions of this software are copyright by UChicago Argonne, LLC and ViennaCL team.

                            -----------------
                  viennaPDE - The Vienna PDE Library
                            -----------------

   Project Head:    Wenyin Wei                   weiwy16@mails.tsinghua.edu.cn

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** \example pde_test.cpp
*
*   In this file, we implemented the concrete numerical schemes taught by Prof. Du Jie and Prof. Yu Hui at Tsinghua University, 2019 autumn.
*  
**/

// System headers
#include <iostream>
#include <fstream>

// ViennaCL headers
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/direct_solve.hpp"
#include "viennacl/linalg/prod.hpp"       //generic matrix-vector product
#include "viennacl/linalg/norm_2.hpp"     //generic l2-norm for vectors
#include "viennacl/linalg/lu.hpp"         //LU substitution routines
#include "viennacl/tools/random.hpp"


// ViennaPDE headers
#include "viennapde/core/varmesh.hpp"
#include "viennapde/core/scheme.hpp"


int main(int argc,char **argv)
{
  


  typedef float   vcl_ScalarT;
  typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
  typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;

  long pixel_index = 0; 
  size_t var_num = 1;
  ssize_t row_num = 1;
  ssize_t column_num = 201;
  vcl_ScalarT TotalTime = 2*3.1415926535;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_varmesh(var_num);

  for (size_t var_index=0; var_index< var_num; var_index++)
  {
    stl_varmesh[var_index].resize(row_num);
    for (ssize_t i=0; i < row_num; ++i) 
      {
        for (ssize_t j= 0; j < column_num; ++j) 
        {
          stl_varmesh[var_index][i].resize(column_num);
          // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
          stl_varmesh[0][i][j] = static_cast<vcl_ScalarT> (std::sin(j*dx)/3*2+1/3);
        }
      }
  }


  viennapde::varmesh<vcl_ScalarT> vcl_varmesh(stl_varmesh);
  viennapde::varmesh<vcl_ScalarT> vcl_varmesh_next(stl_varmesh);

  std::ofstream file;
  for (size_t i = 0; i * dt < TotalTime; i++)
  {
    if (i % 2 ==0) 
    {
      viennapde::scheme::Godunov<vcl_ScalarT>(vcl_varmesh, vcl_varmesh_next, dt, dx);
      viennacl::copy(vcl_varmesh_next.data_[0], stl_varmesh[0]);
    }
    else
    {
      viennapde::scheme::Godunov<vcl_ScalarT>(vcl_varmesh_next, vcl_varmesh, dt, dx);
      viennacl::copy(vcl_varmesh.data_[0], stl_varmesh[0]);
    }
    
    file.open(std::to_string(i) + ".csv");
    if (file.is_open())
      {
        file << stl_varmesh[0][0][0];
        for (size_t j = 1; j < stl_varmesh[0][0].size(); j++)
        {
          
          file << ", "<< stl_varmesh[0][0][j];
        }
        file.close();
      }
  }
  





  /**
  *  That's it.
  **/
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

  return EXIT_SUCCESS;
}

