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
#include "viennapde/test/EXPECT_MESH_EQ.hpp" 

typedef float                         vcl_ScalarT;
typedef viennacl::vector<vcl_ScalarT> vcl_VectorT;
typedef viennacl::matrix<vcl_ScalarT> vcl_MatrixT;
typedef std::vector<std::vector<std::vector<vcl_ScalarT>>>
                                      stl_MeshT;

TEST(MeshClass, copy)
{
  size_t layer_num = 5;
  size_t row_num = 5;
  ssize_t column_num = 10;
  vcl_ScalarT TotalTime = 3.1415926535/10;
  vcl_ScalarT dx = 2*3.1415926535/(column_num-1);
  vcl_ScalarT dt = 0.01;

  stl_MeshT  stl_mesh{layer_num};

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

  // SECTION 02 COPY interface with other classes
  std::vector< viennapde::mesh<vcl_ScalarT>> vcl_mesh;
  vcl_mesh.emplace_back(stl_mesh);
  vcl_mesh.emplace_back(stl_mesh[0].size(), stl_mesh[0][0].size(), stl_mesh.size());
  vcl_mesh.emplace_back(vcl_mesh[0]);
  viennacl::copy(stl_mesh, vcl_mesh[0]);
  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_mesh_new;
  viennacl::copy(vcl_mesh[0], stl_mesh_new);

  EXPECT_MESH_EQ(stl_mesh, stl_mesh_new, "Values diff after copying");

  
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

  std::vector< std::vector< std::vector<vcl_ScalarT> > >  stl_mesh{layer_num};
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::mesh<vcl_ScalarT> 
    vcl_mesh1{stl_mesh}, 
    vcl_mesh2{vcl_mesh1.get_size_num()},
    vcl_mesh3{vcl_mesh1.get_size_num()};
  stl_MeshT stl_mesh_new1, stl_mesh_new2;
  viennacl::copy(vcl_mesh1, stl_mesh_new1);

  vcl_mesh2 = vcl_mesh1 * 2;
  vcl_mesh2 /= 2.0;
  viennacl::copy(vcl_mesh2, stl_mesh_new2);
  EXPECT_MESH_EQ(stl_mesh_new1, stl_mesh_new2, "* & /= failed when *2 /=2.0");


  // The accuracy of + - normally is lower than * -, so users need a higher original value.
  vcl_mesh2 = vcl_mesh1 + 3.0;
  vcl_mesh2 -= 3;
  viennacl::copy(vcl_mesh2, stl_mesh_new2);
  EXPECT_MESH_EQ(stl_mesh_new1, stl_mesh_new2, "+ & -= failed when +3.0 -=3");
  
  // auto vcl_varmat = viennacl::linalg::element_atan(*(vcl_mesh2[0]),*(vcl_mesh1[0]));
  // viennacl::matrix<vcl_ScalarT> vcl_varmat2(row_num, column_num); 
  // vcl_varmat2 = vcl_varmat;
  // viennacl::copy(vcl_varmat2, stl_mesh[0]);
  // for (size_t row_i=0; row_i < row_num; row_i++) {
  //     for (size_t column_i= 0; column_i < column_num; column_i++) {
  //       std::cout << stl_mesh[0][row_i][column_i] << " ";
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

  stl_MeshT  stl_mesh{layer_num}, stl_mesh_new;
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::mesh<vcl_ScalarT> vcl_mesh{stl_mesh};
  vcl_mesh.resize_ptr(2*layer_num);
  vcl_mesh.resize_ptr(layer_num);
  viennacl::copy(vcl_mesh, stl_mesh_new);
  EXPECT_MESH_EQ(stl_mesh, stl_mesh_new, "resize_ptr expand and shrink later, but the content varies.");

  vcl_mesh.clear();
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

  stl_MeshT  stl_mesh{layer_num}, stl_mesh_new;
  for (size_t layer_i=0; layer_i< layer_num; layer_i++) {
    stl_mesh[layer_i].resize(row_num);
    for (size_t row_i=0; row_i < row_num; row_i++) {
      stl_mesh[layer_i][row_i].resize(column_num);
      for (size_t column_i= 0; column_i < column_num; column_i++) {
        // TODO: Hard Code Input, polish the input to make the code can analyze the function automatically.
        stl_mesh[layer_i][row_i][column_i] = static_cast<vcl_ScalarT> (std::sin(column_i*dx)*10+30);
      }
    }
  }

  viennapde::mesh<vcl_ScalarT> vcl_mesh{stl_mesh};
  viennapde::mesh<vcl_ScalarT> vcl_mesh_range{vcl_mesh.begin()+1, vcl_mesh.end()-1};
}



