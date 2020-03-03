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

/** @file viennapde/scheme/pde/WENO.hpp
    @brief Weighted essential non-oscillating methods.
*/
#include <vector>

#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mesh.hpp"
#include "viennapde/core/elem_op.hpp"
#include "viennapde/scheme/pde/scheme_enum.hpp"
#include "viennapde/scheme/pde/flux_func.hpp"
#include "viennapde/scheme/pde/edge_value_estimate.hpp"

// SECTION 01b Declare the image class
namespace viennapde {
namespace scheme {

template <typename NumericT>
/**
 * @brief Impose 1D scheme LaxWendroff + finite volume on the Burger's equation. Requires (0, 2, 0) or more boundary.
 * @param  {const mesh<NumericT> &} iMesh : 
 * @param  {NumericT} dt                  : 
 * @param  {NumericT} dx                  : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> LaxWendroff(const mesh<NumericT> & iMesh, NumericT dt, NumericT dx)
{
    // Calculate \hat{f} at cell edge.
    mesh<NumericT> hatf{iMesh.get_size_num()};
    {
    // Finite Volumn 3 order estimation of the cell edge value.
    auto u_neg = FV_EdgeEstimate<NumericT, EdgeFaceOrien::NEG>(iMesh);
    auto u_pos = FV_EdgeEstimate<NumericT, EdgeFaceOrien::POS>(iMesh);

    hatf = FluxFunction<NumericT, FluxType::LaxWendroff_flux>(u_neg, u_pos);
    }

    // Calculate \hat{f}_{j+1/2} - \hat{f}_{j-1/2}
    mesh<NumericT> hatfdiff{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 0) = -1; tKernel(0, 1) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,0}, {0,1}};
    hatfdiff = convolve(hatf, tKernel, tROIrc_vec);
    }
    // Time Forward, attention this is a temporarlly extended mesh to suffice the need of periodic B.C..
    return iMesh - hatfdiff * (dt / dx);
}

} //namespace viennapde::scheme
} //namespace viennapde
