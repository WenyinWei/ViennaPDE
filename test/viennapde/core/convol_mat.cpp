/* =========================================================================
   Copyright (c) 2010-2020, Department of Engineering Physics,
                            Tsinghua University, Beijing, China.

   Portions of this software are copyright by UChicago Argonne, LLC and ViennaCL team.

                            -----------------
                  viennapde - The Vienna PDE Library
                            -----------------

   Project Head:    Wenyin Wei                   weiwy16@mails.tsinghua.edu.cn

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file test/viennapde/core/convol.hpp
    @brief TEST: Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library. The 2D-style data handling utilizes the potential of GPU to its extent.
*/

#include <cmath>

#include "gtest/gtest.h"
#include "viennapde/test/EXPECT_MESH_EQ.hpp"

#include "viennapde/core/convol_mat.hpp"

typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef viennapde::mesh<vcl_ScalarT>  vcl_MeshT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;

class ConvolMatonMat : public ::testing::Test {
 protected:
  void SetUp() override { 
    for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
      stl_mesh[layer_i].resize(row_num);
      for (size_t row_i=0; row_i < row_num; row_i++) {
        stl_mesh[layer_i][row_i].resize(column_num);
        for (size_t column_i= 0; column_i < column_num; column_i++) {
          // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
          stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1.0/3);
        }
      }
    }

    vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size()); //vcl_mesh_vec[0]
    vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size()); //vcl_mesh_vec[1]
    viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
    
    vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[0]
    vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[1]
    for (size_t row_i=0; row_i < kernel_size1; row_i++) {
      for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
        vcl_ker_vec[0](row_i, column_i) = std::sin(column_i*dx)/3*2+1/3;
        vcl_ker_vec[1](row_i, column_i) = - vcl_ker_vec[0](row_i, column_i);
      }
    }
  }

  // void TearDown() override {}

  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 1;

  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};
  std::vector< vcl_MeshT > vcl_mesh_vec;
  std::vector< vcl_MatrixT> vcl_ker_vec;
  std::vector< stl_MeshT > stl_mesh_vec_new{2};
};


#define CONVOL_MATONMAT(CONVOLTYPE)          \
TEST_F(ConvolMatonMat, CONVOLTYPE)          \
{                                           \
  constexpr auto convolT = viennapde::ConvolutionType::CONVOLTYPE;                    \
                                                                                      \
  viennapde::convolve<vcl_ScalarT, convolT>(                                          \
      *(vcl_mesh_vec[0][0]), vcl_ker_vec[0], *(vcl_mesh_vec[1][0]));                  \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                               \
  viennapde::convolve<vcl_ScalarT, convolT>(                                          \
      *(vcl_mesh_vec[0][0]), vcl_ker_vec[1], *(vcl_mesh_vec[1][0]));                  \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                               \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                   \
    "Value results are not right by a kernel and its negative corresponding one.");   \
                                                                                      \
  *(vcl_mesh_vec[1][0]) = viennapde::convolve<vcl_ScalarT, convolT>(                  \
      *(vcl_mesh_vec[0][0]), vcl_ker_vec[0]);                                         \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                               \
  *(vcl_mesh_vec[1][0]) = viennapde::convolve<vcl_ScalarT, convolT>(                  \
      *(vcl_mesh_vec[0][0]), vcl_ker_vec[1]);                                         \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                               \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                   \
    "Value results are not right by a kernel and its negative corresponding one.");   \
}

CONVOL_MATONMAT(INNER);
CONVOL_MATONMAT(EQUIV);
CONVOL_MATONMAT(OUTER);

#undef CONVOL_MATONMAT
