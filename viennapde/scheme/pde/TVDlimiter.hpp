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

/** @file viennapde/scheme/pde/TVDlimiter.hpp
    @brief Total variation diminishing limiter.
*/

#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mesh.hpp"
#include "viennapde/core/elem_op.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

// TODO Minimod function should be more flexible, e.g., support number-casual input as an arg list.
template <typename NumericT>
mesh<NumericT> minimod(const mesh<NumericT> & iMesh1, const mesh<NumericT> & iMesh2, const mesh<NumericT> & iMesh3)
{
    assert(iMesh1.get_size_num()==iMesh2.get_size_num() && iMesh2.get_size_num()==iMesh3.get_size_num());
    return elem_geq0(iMesh1) * elem_geq0(iMesh2) * elem_geq0(iMesh3) * elem_min( elem_min(iMesh1, iMesh2), iMesh3)
        +  elem_leq0(iMesh1) * elem_leq0(iMesh2) * elem_leq0(iMesh3) * elem_max( elem_max(iMesh1, iMesh2), iMesh3);
}

} //namespace viennapde::scheme
} //namespace viennapde
