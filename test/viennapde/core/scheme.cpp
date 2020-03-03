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
typedef viennapde::mesh<vcl_ScalarT>  vcl_MeshT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;

class Scheme1D : public ::testing::Test {
 protected:
  void SetUp() override { 
    for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
      stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
    for (size_t column_i= 0; column_i < column_num; column_i++) {
      stl_mesh[layer_i][row_i][column_i] =  1.0/3 + std::sin(column_i*dx)*2/3;
    }}}
  viennacl::copy(stl_mesh, vcl_mesh);
  }

  void TearDown() override {
    const size_t by = mb.get_marginy();
    for (size_t i = 0; i * dt < TotalTime; i++) {
      mb.refresh();
      viennacl::copy(vcl_mesh, stl_mesh);
      file.open(data_path + std::to_string(i) + ".csv");
      if (file.is_open())
      {
        file << stl_mesh[0][0][by];
        for (size_t column_i = by+1; column_i < stl_mesh[0][0].size()-by; column_i++)
          file << ", "<< stl_mesh[0][0][column_i];
        file << "\n";
        file.close();
      } else { std::cerr << "Fail to open files to store data.\n"; };
      vcl_mesh = scheme(vcl_mesh, dt, dx);
    }
  }

  size_t layer_num = 1;
  size_t row_num = 1;
  size_t column_num = 60;
  vcl_ScalarT TotalTime = 5 * M_PI;
  vcl_ScalarT dx = 2* M_PI / column_num;
  vcl_ScalarT dt = 0.02;

  stl_MeshT stl_mesh{layer_num}, stl_varmesh_new;
  vcl_MeshT  vcl_mesh{row_num, column_num, layer_num};
  viennapde::meshb<vcl_ScalarT> mb{vcl_mesh};

  std::function<vcl_MeshT (const vcl_MeshT & , vcl_ScalarT, vcl_ScalarT )> scheme;

  std::ofstream file; 
  std::string data_path;
};

TEST_F(Scheme1D, Warmup)     
{                                       
  mb.extend(0, 1, 0);
  data_path = DATAPATH_SCHEME1D_WARMUP ;

  scheme = viennapde::scheme::Godunov<vcl_ScalarT>;
} 

TEST_F(Scheme1D, Godunov)     
{                                       
  mb.extend(0, 1, 0);
  data_path = DATAPATH_SCHEME1D_GODUNOV ;

  scheme = viennapde::scheme::Godunov<vcl_ScalarT>;
} 

TEST_F(Scheme1D, Godunov_RK3order)     
{                                       
  mb.extend(0, 3, 0);
  
  data_path = DATAPATH_SCHEME1D_GODUNOV_RK3ORDER ;
  
  using namespace std::placeholders;  // for _1, _2, _3...
  scheme = std::bind(viennapde::scheme::RK3order<vcl_ScalarT>, _1, _2, _3, viennapde::scheme::Godunov<vcl_ScalarT>);
} 

TEST_F(Scheme1D, LaxWendroff)     
{                                       
  mb.extend(0, 2, 0);
  data_path = DATAPATH_SCHEME1D_LAXWENDROFF ;

  scheme = viennapde::scheme::LaxWendroff<vcl_ScalarT>;
} 

TEST_F(Scheme1D, WENO_LaxFriedrichs)     
{                                       
  mb.extend(0, 3, 0);
  data_path = DATAPATH_SCHEME1D_WENO_LAXFRIEDRICHS ;

  scheme = viennapde::scheme::WENO<vcl_ScalarT, viennapde::FluxType::LaxFriedrichs>;
} 

