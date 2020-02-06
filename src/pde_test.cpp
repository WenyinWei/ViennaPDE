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
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/scheme.hpp"


int main(int argc,char **argv)
{
  typedef float                         vcl_ScalarT;
  typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
  typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;

  long pixel_index = 0; 
  size_t layer_num = 5;
  size_t row_num = 50;
  ssize_t column_num = 50;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_varmesh(layer_num);

  for (size_t layer_i=0; layer_i< layer_num; layer_i++)
  {
    stl_varmesh[layer_i].resize(row_num);
    for (ssize_t row_i=0; row_i < row_num; ++row_i) 
      {
        stl_varmesh[layer_i][row_i].resize(column_num);
        for (ssize_t column_i= 0; column_i < column_num; ++column_i) 
        {
          // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
          stl_varmesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
        }
      }
  }


  viennapde::Varmesh<vcl_ScalarT> vcl_varmesh{stl_varmesh};
  viennapde::Varmesh<vcl_ScalarT> vcl_varmesh_next{stl_varmesh};

  std::ofstream file;
  for (size_t row_i = 0; row_i * dt < TotalTime; row_i++)
  {
    if (row_i % 2 ==0) 
    {
      viennapde::scheme::Godunov<vcl_ScalarT>(vcl_varmesh, vcl_varmesh_next, dt, dx);
      viennacl::copy(vcl_varmesh_next, stl_varmesh);
    }
    else
    {
      viennapde::scheme::Godunov<vcl_ScalarT>(vcl_varmesh_next, vcl_varmesh, dt, dx);
      viennacl::copy(vcl_varmesh, stl_varmesh);
    }
    
    file.open(std::to_string(row_i) + ".csv");
    if (file.is_open())
      {
        file << stl_varmesh[0][0][0];
        for (size_t column_i = 1; column_i < stl_varmesh[0][0].size(); column_i++)
        {
          
          file << ", "<< stl_varmesh[0][0][column_i];
        }
        file.close();
      }
  }
  viennapde::scheme::RK3order<vcl_ScalarT>(vcl_varmesh, vcl_varmesh_next, dt, dx, viennapde::scheme::Godunov<vcl_ScalarT>);
  





  std::cout << "I am fine 5\n";
  /**
  *  That's it.
  **/
  std::cout << "!!!! TUTORIAL COMPLETED SUCCESSFULLY !!!!" << std::endl;

  return EXIT_SUCCESS;
}

