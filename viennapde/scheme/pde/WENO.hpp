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
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mesh.hpp"
#include "viennapde/core/elem_op.hpp"
#include "viennapde/scheme/pde/scheme_enum.hpp"
#include "viennapde/scheme/pde/flux_func.hpp"
#include "viennapde/scheme/pde/edge_value_estimate.hpp"
#include "viennapde/scheme/pde/TVDlimiter.hpp"

// SECTION 01b Declare the image class
namespace viennapde {
namespace scheme {

template <typename NumericT, FluxType fluxT>
/**
 * @brief Impose 1D scheme LaxFriedrichs + finite volume on the Burger's equation. Requires (0, 3, 0) or more boundary.
 * @param  {const mesh<NumericT> &} iMesh : 
 * @param  {NumericT} dt                  : 
 * @param  {NumericT} dx                  : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> WENO(const mesh<NumericT> & iMesh, NumericT dt, NumericT dx)
{
    // Calculate u_{j+1/2}^-, u_{j+1/2}^+ for f(u) estimation at cell edge.
    mesh<NumericT> hatf{iMesh.get_size_num()};
    {   
    // Finite Volumn 5 order estimation of the cell edge value.
    auto u_neg = WENO_g<NumericT, EdgeFaceOrien::NEG>(iMesh);
    auto u_pos = WENO_g<NumericT, EdgeFaceOrien::POS>(iMesh);
    
    {// Minimod function as TVD limiter
    viennacl::matrix<NumericT> tKernel3{1, 3}; tKernel3(0,1)=-1.0; tKernel3(0,2)=1.0; 
    std::vector<cord2<GridIntT>> tROIrc_vec3{{0,1}, {0,2}};
    auto u_jp1_m_u_j = convolve(iMesh, tKernel3, tROIrc_vec3);
    
    viennacl::matrix<NumericT> tKernel4{1, 3}; tKernel4(0,1)= 1.0; tKernel4(0,0)=-1.0; 
    std::vector<cord2<GridIntT>> tROIrc_vec4{{0,0}, {0,1}};
    auto u_j_m_u_jm1 = convolve(iMesh, tKernel4, tROIrc_vec4);
    
    viennacl::matrix<NumericT> tKernel5{1, 3}; tKernel5(0,0)=1.0;  
    std::vector<cord2<GridIntT>> tROIrc_vec5{{0,0}};
    auto u_jm1_pos   = convolve(u_pos, tKernel5, tROIrc_vec5);
    
    // auto tilde_u_mod = minimod(u_neg - iMesh,       u_jp1_m_u_j, u_j_m_u_jm1);
    // auto ttilde_u_mod= minimod(iMesh - u_jm1_pos,   u_jp1_m_u_j, u_j_m_u_jm1);
    
    viennacl::matrix<NumericT> tKernel6{1, 3}; tKernel6(0,2)=1.0;  
    std::vector<cord2<GridIntT>> tROIrc_vec6{{0,2}};
    // u_neg => u_neg(mod)  u_pos => u_pos(mod)
    u_neg = iMesh          + minimod(u_neg - iMesh,     u_jp1_m_u_j, u_j_m_u_jm1);
    u_pos = convolve(iMesh - minimod(iMesh - u_jm1_pos, u_jp1_m_u_j, u_j_m_u_jm1), tKernel6, tROIrc_vec6);
    }

    hatf = FluxFunction<NumericT, fluxT>(u_neg, u_pos);
    }

    // Calculate \hat{f}_{j+1/2} - \hat{f}_{j-1/2}
    mesh<NumericT> hatfdiff{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 0) = -1; tKernel(0, 1) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,0}, {0,1}};
    hatfdiff = convolve(hatf, tKernel, tROIrc_vec);
    }
    // Time Forward
    return iMesh - hatfdiff * (dt / dx);
}


} //namespace viennapde::scheme
} //namespace viennapde
