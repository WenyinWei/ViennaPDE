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

/** @file viennapde/scheme/RK3order.hpp
    @brief Implementation of Runge-Kuta 3 order scheme.
*/

#include <functional>

#include "viennapde/core/mesh.hpp"
#include "viennapde/core/meshb.hpp"
#include "viennapde/core/convol_mesh.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

template <typename NumericT>
mesh<NumericT> RK3order(const mesh<NumericT> & iMesh, NumericT dt, NumericT dx, 
    std::function<mesh<NumericT> (const mesh<NumericT> & , NumericT, NumericT )> scheme)
{
    std::vector< mesh<NumericT> > tMesh;
    for (size_t i = 0; i < 2; i++)
        tMesh.emplace_back(iMesh.get_size_num());
    tMesh[0] = scheme(iMesh,    dt, dx);
    tMesh[1] = scheme(tMesh[0], dt, dx);
    tMesh[1] = iMesh * (3.0/4) + tMesh[1] * (1.0/4);
    return iMesh * (1.0/3) + scheme(tMesh[1], dt, dx) * (2.0/3);
}

} //namespace viennapde::scheme
} //namespace viennapde
