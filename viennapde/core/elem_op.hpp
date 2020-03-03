#pragma once
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

/** @file viennapde/scheme/elem_op.hpp
    @brief Implementation of concrete scheme for PDE.
    TODO There may be abundant includes, remove them.
*/
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>
#include <iostream>

#include "viennacl/linalg/matrix_operations.hpp"
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/meshb.hpp"
#include "viennapde/core/convol_mesh.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
template <typename NumericT>
mesh<NumericT> elem_geq0(const mesh<NumericT> & iMesh)
{
    mesh<NumericT> tMesh{iMesh.get_size_num()};
    viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(iMesh.get_row_num(), iMesh.get_column_num(), 0.5);
    for (size_t i = 0; i < iMesh.get_layer_num(); i++)
        *(tMesh[i]) = M_1_PI * viennacl::linalg::element_atan(*(iMesh[i]) * 100000.0) + scalar_mat;
    return tMesh;
}

template <typename NumericT>
mesh<NumericT> elem_leq0(const mesh<NumericT> & iMesh)
{
    mesh<NumericT> tMesh{iMesh.get_size_num()};
    viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(iMesh.get_row_num(), iMesh.get_column_num(), 0.5);
    for (size_t i = 0; i < iMesh.get_layer_num(); i++)
        *(tMesh[i]) = M_1_PI * viennacl::linalg::element_atan(*(iMesh[i]) * -100000.0) + scalar_mat;
    return tMesh;
}

template <typename NumericT>
mesh<NumericT> elem_abs(const mesh<NumericT> & iMesh)
{
    mesh<NumericT> tMesh{iMesh.get_size_num()};
    for (size_t i = 0; i < iMesh.get_layer_num(); i++)
        *(tMesh[i]) = viennacl::linalg::element_prod( 
                        M_2_PI * viennacl::linalg::element_atan(*(iMesh[i]) * 100000.0), *(iMesh[i]));
    return tMesh;
}

template <typename NumericT>
mesh<NumericT> elem_max(const mesh<NumericT> & iMeshlhs, const mesh<NumericT> & iMeshrhs)
{
    assert(iMeshlhs.get_size_num() == iMeshrhs.get_size_num());
    return iMeshlhs + (iMeshrhs-iMeshlhs) * elem_geq0(iMeshrhs-iMeshlhs);
}

template <typename NumericT>
mesh<NumericT> elem_min(const mesh<NumericT> & iMeshlhs, const mesh<NumericT> & iMeshrhs)
{
    assert(iMeshlhs.get_size_num() == iMeshrhs.get_size_num());
    return iMeshlhs + (iMeshrhs-iMeshlhs) * elem_leq0(iMeshrhs-iMeshlhs);
}
} //namespace viennapde
