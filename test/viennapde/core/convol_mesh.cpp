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

/** @file test/viennapde/core/convol_mesh.cpp
    @brief TEST: Convolve operation on mesh class. All functions in this file are overloaded of viennapde::convolve working on the mesh class. 
*/

#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "viennapde/test/EXPECT_MESH_EQ.hpp"

#include "viennapde/core/convol_mesh.hpp"


typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef viennapde::mesh<vcl_ScalarT>  vcl_MeshT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;

class ConvolMatonMesh : public ::testing::Test {
 protected:
  void SetUp() override { 
    for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
      stl_mesh[layer_i].resize(row_num);
      for (size_t row_i=0; row_i < row_num; row_i++) {
        stl_mesh[layer_i][row_i].resize(column_num);
        for (size_t column_i= 0; column_i < column_num; column_i++) {
          // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
          stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
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
  size_t layer_num = 3;

  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};
  std::vector< vcl_MeshT > vcl_mesh_vec;
  std::vector< vcl_MatrixT> vcl_ker_vec;
  std::vector< stl_MeshT > stl_mesh_vec_new{2};
};

#define CONVOLMATONMESH(CONVOLTYPE)     \
TEST_F(ConvolMatonMesh, CONVOLTYPE)     \
{                                       \
  constexpr auto convolT = viennapde::ConvolutionType::CONVOLTYPE;                  \
                                                                                    \
  viennapde::convolve<vcl_ScalarT, convolT>(                                        \
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);                            \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                             \
  viennapde::convolve<vcl_ScalarT, convolT>(                                        \
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);                            \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                             \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                 \
    "Value results are not right by a kernel and its negative corresponding one."); \
                                                                                    \
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(                      \
      vcl_mesh_vec[0], vcl_ker_vec[0]);                                             \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                             \
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(                      \
      vcl_mesh_vec[0], vcl_ker_vec[1]);                                             \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                             \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                 \
    "Value results are not right by a kernel and its negative corresponding one."); \
} 

CONVOLMATONMESH(INNER);
CONVOLMATONMESH(EQUIV);
CONVOLMATONMESH(OUTER);

#undef CONVOLMATONMESH


class ConvolMeshonMesh : public ::testing::Test {
 protected:
  void SetUp() override { 
    for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
      stl_mesh[layer_i].resize(row_num);
      for (size_t row_i=0; row_i < row_num; row_i++) {
        stl_mesh[layer_i][row_i].resize(column_num);
        for (size_t column_i= 0; column_i < column_num; column_i++) {
          // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
          stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
        }
      }
    }

    vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size()); //vcl_mesh_vec[0]
    vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size()); //vcl_mesh_vec[1]
    viennacl::copy(stl_mesh, vcl_mesh_vec[0]);

    vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[0]
    vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[1]
    for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)
      for (size_t row_i=0; row_i < kernel_size1; row_i++) {
          for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
              (*(vcl_ker_vec[0][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
              (*(vcl_ker_vec[1][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
          }
      }
  }

  // void TearDown() override {}

  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3, kernel_size3 = 3;

  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};
  std::vector< vcl_MeshT > vcl_mesh_vec;
  std::vector< vcl_MeshT > vcl_ker_vec;
  std::vector< stl_MeshT > stl_mesh_vec_new{2};
};

#define CONVOLMESHONMESH(CONVOLTYPE)    \
TEST_F(ConvolMeshonMesh, CONVOLTYPE)    \
{                                       \
  constexpr auto convolT = viennapde::ConvolutionType::CONVOLTYPE;                  \
                                                                                    \
  viennapde::convolve<vcl_ScalarT, convolT>(                                        \
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);                            \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                             \
  viennapde::convolve<vcl_ScalarT, convolT>(                                        \
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);                            \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                             \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                 \
    "Value results are not right by a kernel and its negative corresponding one."); \
                                                                                    \
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(                      \
      vcl_mesh_vec[0], vcl_ker_vec[0]);                                             \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[0]);                             \
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(                      \
      vcl_mesh_vec[0], vcl_ker_vec[1]);                                             \
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_new[1]);                             \
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_new[0], stl_mesh_vec_new[1],                 \
    "Value results are not right by a kernel and its negative corresponding one."); \
                                                                                    \
  vcl_ScalarT result = 0;                                                           \
  for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)                       \
    for (size_t row_i=0; row_i < kernel_size1; row_i++)                             \
        for (size_t column_i= 0; column_i < kernel_size2; column_i++) {             \
            result += (*(vcl_mesh_vec[0][layer_i]))(row_i, column_num) * (*(vcl_ker_vec[0][layer_i]))(row_i, column_i);\
        }                                                                           \
  if constexpr (convolT == viennapde::ConvolutionType::INNER)                       \
  {                                                                                 \
    if constexpr (std::is_same<vcl_ScalarT, float>::value) {                        \
      EXPECT_FLOAT_EQ(result, stl_mesh_vec_new[0][0][0][0]) << "The result is not right under a 3*3*3 mesh convolved by a 3*3*3 mesh ker.\n";\
    } else if constexpr (std::is_same<vcl_ScalarT, float>::value) {                 \
      EXPECT_DOUBLE_EQ(result, stl_mesh_vec_new[0][0][0][0]) << "The result is not right under a 3*3*3 mesh convolved by a 3*3*3 mesh ker.\n";\
    } else {FAIL() << "Specify your wanted scalar value type to test, float or double?";}\
  }                                                                                 \
} 

CONVOLMESHONMESH(INNER);
CONVOLMESHONMESH(EQUIV);
CONVOLMESHONMESH(OUTER);

#undef CONVOLMESHONMESH