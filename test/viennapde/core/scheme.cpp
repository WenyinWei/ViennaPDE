/* =========================================================================
   Copyright (c) 2016-2020, Department of Engineering Physics,
                            Tsinghua University, Beijing, China.

   Portions of this software are copyright by UChicago Argonne, LLC and ViennaCL team.

                            -----------------
                  ViennaPDE - The Vienna PDE Library
                            -----------------

   Project Head:    Wenyin Wei                   weiwy16@mails.tsinghua.edu.cn

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file test/viennapde/core/scheme.hpp
    @brief TEST: Implementation of concrete scheme for PDE.
*/

#include <cmath>
#include <functional>
#include <filesystem>

#include "gtest/gtest.h"

#include "viennapde/core/scheme.hpp"

typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;
std::string data_folder = "./data/scheme/Godunov/";

TEST(Scheme, Godunov)
{
  size_t layer_num = 1;
  size_t row_num = 1;
  size_t column_num = 100;
  vcl_ScalarT TotalTime = 2 * M_PI;
  vcl_ScalarT dx = 2* M_PI / column_num;
  vcl_ScalarT dt = 0.02;

  stl_MeshT stl_varmesh{layer_num}, stl_varmesh_after;

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_varmesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_varmesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_varmesh[layer_i][row_i][column_i] =  2.0/3 + std::sin(column_i*dx)*1/3;
      }
    }
  }

  viennapde::Varmesh<vcl_ScalarT>  vcl_varmesh{row_num, column_num, layer_num};
  viennacl::copy(stl_varmesh, vcl_varmesh);

  viennapde::meshb mb(vcl_varmesh);
  mb.extend(0, 1, 0);
  mb.refreshBC();
  
  
  std::ofstream file; 
  const size_t by = mb.get_marginy();
  for (size_t i = 0; i * dt < TotalTime; i++)
  {
    mb.refreshBC();
    viennacl::copy(vcl_varmesh, stl_varmesh);
    file.open(data_folder + std::to_string(i) + ".csv");
    if (file.is_open())
    {
      file << stl_varmesh[0][0][by];
      for (size_t column_i = by+1; column_i < stl_varmesh[0][0].size()-by; column_i++)
        file << ", "<< stl_varmesh[0][0][column_i];
      file.close();
    } else { std::cerr << "Fail to open files to store data.\n"; };
    vcl_varmesh = viennapde::scheme::Godunov(vcl_varmesh, dt, dx);
    // vcl_varmesh = viennapde::scheme::RK3order(vcl_varmesh, dt, dx,
                                                // viennapde::scheme::Godunov<vcl_ScalarT>);
  }

//   EXPECT_MESH_EQ(stl_varmesh, stl_varmesh_after, "Values diff after copying");
}
