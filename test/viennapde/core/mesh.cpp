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

/** @file test/viennapde/core/mesh.cpp
    @brief TEST: Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
*/
#include <cmath>

#include "gtest/gtest.h"

#include "viennapde/core/mesh.hpp"

typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;

void EXPECT_MESH_EQ(stl_MeshT & A, stl_MeshT & B, std::string err_info) {
  for (size_t layer_i=0; layer_i< A.size(); layer_i++) {
    for (size_t row_i=0; row_i < A[0].size(); row_i++) {
      for (size_t column_i= 0; column_i < A[0][0].size(); column_i++) {
        if constexpr (std::is_same<vcl_ScalarT, float>::value) {
          EXPECT_FLOAT_EQ( A[layer_i][row_i][column_i], B[layer_i][row_i][column_i]) << err_info;
        } else if constexpr (std::is_same<vcl_ScalarT, double>::value) {
          EXPECT_DOUBLE_EQ(A[layer_i][row_i][column_i], B[layer_i][row_i][column_i]) << err_info;
        } else {FAIL() << "Specify your wanted scalar value type to test, float or double?";}
      }
    }
  }
}

TEST(MeshClass, copy)
{
  size_t layer_num = 5;
  size_t row_num = 5;
  ssize_t column_num = 10;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  stl_MeshT  stl_varmesh{layer_num};

  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_varmesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_varmesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_varmesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)/3*2+1/3);
      }
    }
  }

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::Varmesh<vcl_ScalarT>> vcl_varmesh;
  vcl_varmesh.emplace_back(stl_varmesh);
  vcl_varmesh.emplace_back(stl_varmesh[0].size(), stl_varmesh[0][0].size(), stl_varmesh.size());
  vcl_varmesh.emplace_back(vcl_varmesh[0]);
  viennacl::copy(stl_varmesh, vcl_varmesh[0]);
  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_varmesh_after;
  viennacl::copy(vcl_varmesh[0], stl_varmesh_after);

  EXPECT_MESH_EQ(stl_varmesh, stl_varmesh_after, "Values diff after copying");

  
}


/*===== SECTION Operator Overloading ======================================== */

TEST(MeshClass, operators)
{
  size_t layer_num = 5;
  size_t row_num = 5;
  ssize_t column_num = 10;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_varmesh{layer_num};
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_varmesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_varmesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_varmesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::Varmesh<vcl_ScalarT> 
    vcl_varmesh1{stl_varmesh}, 
    vcl_varmesh2{vcl_varmesh1.get_size_num()},
    vcl_varmesh3{vcl_varmesh1.get_size_num()};
  stl_MeshT stl_varmesh_new1, stl_varmesh_new2;
  viennacl::copy(vcl_varmesh1, stl_varmesh_new1);

  vcl_varmesh2 = vcl_varmesh1 * 2;
  vcl_varmesh2 /= 2.0;
  viennacl::copy(vcl_varmesh2, stl_varmesh_new2);
  EXPECT_MESH_EQ(stl_varmesh_new1, stl_varmesh_new2, "* & /= failed when *2 /=2.0");


  // The accuracy of + - normally is lower than * -, so users need a higher original value.
  vcl_varmesh2 = vcl_varmesh1 + 3.0;
  vcl_varmesh2 -= 3;
  viennacl::copy(vcl_varmesh2, stl_varmesh_new2);
  EXPECT_MESH_EQ(stl_varmesh_new1, stl_varmesh_new2, "+ & -= failed when +3.0 -=3");
  
  // auto vcl_varmat = viennacl::linalg::element_atan(*(vcl_varmesh2[0]),*(vcl_varmesh1[0]));
  // viennacl::matrix<vcl_ScalarT> vcl_varmat2(row_num, column_num); 
  // vcl_varmat2 = vcl_varmat;
  // viennacl::copy(vcl_varmat2, stl_varmesh[0]);
  // for (size_t row_i=0; row_i < row_num; row_i++) {
  //     for (size_t column_i= 0; column_i < column_num; column_i++) {
  //       std::cout << stl_varmesh[0][row_i][column_i] << " ";
  //     }
  //     std::cout << "\n";
  //   }
}




/*===== SECTION Operator Overloading ======================================== */

TEST(MeshClass, resize_ptr_and_clear)
{
  size_t layer_num = 5;
  size_t row_num = 5;
  ssize_t column_num = 10;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  stl_MeshT  stl_varmesh{layer_num}, stl_varmesh_new;
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_varmesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_varmesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_varmesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::Varmesh<vcl_ScalarT> vcl_varmesh{stl_varmesh};
  vcl_varmesh.resize_ptr(2*layer_num);
  vcl_varmesh.resize_ptr(layer_num);
  viennacl::copy(vcl_varmesh, stl_varmesh_new);
  EXPECT_MESH_EQ(stl_varmesh, stl_varmesh_new, "resize_ptr expand and shrink later, but the content varies.");

  vcl_varmesh.clear();
}


/*===== SECTION Operator Overloading ======================================== */

TEST(MeshClass, range_CTOR)
{
  size_t layer_num = 5;
  size_t row_num = 5;
  ssize_t column_num = 10;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  stl_MeshT  stl_varmesh{layer_num}, stl_varmesh_new;
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_varmesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_varmesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_varmesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::Varmesh<vcl_ScalarT> vcl_varmesh{stl_varmesh};
  viennapde::Varmesh<vcl_ScalarT> vcl_varmesh_range{vcl_varmesh.begin()+1, vcl_varmesh.end()-1};
}



