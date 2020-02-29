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
    TODO Not yet and verified.
*/



#include <cmath>
#include <vector>

#include "gtest/gtest.h"
#include "viennapde/test/EXPECT_MESH_EQ.hpp"

#include "viennapde/core/convol_mesh.hpp"



typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;


TEST(Convol, MatonMesh_INNER)
{
  constexpr auto convolT = viennapde::ConvolutionType::INNER;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennacl::matrix<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[1]
  for (size_t row_i=0; row_i < kernel_size1; row_i++) {
    for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
      vcl_ker_vec[0](row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      vcl_ker_vec[1](row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
    }
  }

  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
}


TEST(Convol, MatonMesh_OUTER)
{
  constexpr auto convolT = viennapde::ConvolutionType::OUTER;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennacl::matrix<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[1]
  for (size_t row_i=0; row_i < kernel_size1; row_i++) {
    for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
      vcl_ker_vec[0](row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      vcl_ker_vec[1](row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
    }
  }

  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
}


TEST(Convol, MatonMesh_EQUIV)
{
  constexpr auto convolT = viennapde::ConvolutionType::EQUIV;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennacl::matrix<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2); // vcl_ker_vec[1]
  for (size_t row_i=0; row_i < kernel_size1; row_i++) {
    for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
      vcl_ker_vec[0](row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      vcl_ker_vec[1](row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
    }
  }

  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
}



TEST(Convol, MeshonMesh_OUTER)
{
  constexpr auto convolT = viennapde::ConvolutionType::OUTER;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3, kernel_size3 = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennapde::mesh<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[1]
  for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)
    for (size_t row_i=0; row_i < kernel_size1; row_i++) {
        for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
            (*(vcl_ker_vec[0][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
            (*(vcl_ker_vec[1][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
        }
    }
  


  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
}

TEST(Convol, MeshonMesh_EQUIV)
{
  constexpr auto convolT = viennapde::ConvolutionType::EQUIV;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3, kernel_size3 = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennapde::mesh<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[1]
  for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)
    for (size_t row_i=0; row_i < kernel_size1; row_i++) {
        for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
            (*(vcl_ker_vec[0][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
            (*(vcl_ker_vec[1][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
        }
    }
  


  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
}

TEST(Convol, MeshonMesh_INNER)
{
  constexpr auto convolT = viennapde::ConvolutionType::INNER;
  size_t row_num = 3, kernel_size1 = 3;
  size_t column_num = 3, kernel_size2 = 3;
  size_t layer_num = 3, kernel_size3 = 3;
  vcl_ScalarT dx = 1.0;
  
  stl_MeshT  stl_mesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }
  

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh_vec;
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh_vec.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  viennacl::copy(stl_mesh, vcl_mesh_vec[0]);
  
  std::vector<viennapde::mesh<vcl_ScalarT>> vcl_ker_vec;
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[0]
  vcl_ker_vec.emplace_back(kernel_size1, kernel_size2, kernel_size3); // vcl_ker_vec[1]
  for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)
    for (size_t row_i=0; row_i < kernel_size1; row_i++) 
        for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
            (*(vcl_ker_vec[0][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
            (*(vcl_ker_vec[1][layer_i]))(row_i, column_i) = static_cast<vcl_ScalarT> (-std::sin(column_i*dx)/3*2-1/3);
        }

  std::vector< stl_MeshT >  stl_mesh_vec_after{2};
  
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1], vcl_mesh_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");
  
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[0]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[0]);
  vcl_mesh_vec[1] = viennapde::convolve<vcl_ScalarT, convolT>(
      vcl_mesh_vec[0], vcl_ker_vec[1]);
  viennacl::copy(vcl_mesh_vec[1], stl_mesh_vec_after[1]);
  EXPECT_MESH_NEGATIVE_EQ(stl_mesh_vec_after[0], stl_mesh_vec_after[1], 
    "Value results are not right by a kernel and its negative corresponding one.");

  float result = 0;
  for (size_t layer_i = 0; layer_i < kernel_size3; layer_i++)
    for (size_t row_i=0; row_i < kernel_size1; row_i++) 
        for (size_t column_i= 0; column_i < kernel_size2; column_i++) {
            result += (*(vcl_mesh_vec[0][layer_i]))(row_i, column_num) * (*(vcl_ker_vec[0][layer_i]))(row_i, column_i);
        }
  EXPECT_EQ(result, stl_mesh_vec_after[0][0][0][0]) << "The result is not right under a 3*3*3 mesh convolved by a 3*3*3 mesh ker.\n";
}