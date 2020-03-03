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

/** @file viennapde/core/scheme.hpp
    @brief Implementation of concrete scheme for PDE.
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
enum FluxType {
    LaxFriedrichs,
    LaxWendroff,
    Godunov
};

// TODO: Change it to operator overloading
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

template <typename NumericT>
mesh<NumericT> minimod(const mesh<NumericT> & iMesh1, const mesh<NumericT> & iMesh2, const mesh<NumericT> & iMesh3)
{
    assert(iMesh1.get_size_num()==iMesh2.get_size_num() && iMesh2.get_size_num()==iMesh3.get_size_num());
    return elem_geq0(iMesh1) * elem_geq0(iMesh2) * elem_geq0(iMesh3) * elem_min( elem_min(iMesh1, iMesh2), iMesh3)
        +  elem_leq0(iMesh1) * elem_leq0(iMesh2) * elem_leq0(iMesh3) * elem_max( elem_max(iMesh1, iMesh2), iMesh3);
}

namespace scheme
{
// This func requires extend_boundary(0,1,1);
template <typename NumericT>
/**
 * @brief Impose 1D scheme Godunov based on the shock and rarefaction analysis on the Burger's equation. Requires (0, 1, 0) or more boundary.
 * @param  {const mesh<NumericT> &} iMesh : 
 * @param  {NumericT} dt                  : 
 * @param  {NumericT} dx                  : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> Godunov(const mesh<NumericT> & iMesh, NumericT dt, NumericT dx)
{
    // Calculate f(u^*)
    mesh<NumericT> hatf{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 2) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,2}};
    mesh<NumericT>         ur = convolve<NumericT, ConvolutionType::EQUIV>(iMesh, tKernel, tROIrc_vec);
    mesh<NumericT> const & ul = iMesh;
    mesh<NumericT> ustar 
        = (ur + (ul-ur) * elem_geq0(ul+ur)) * ((NumericT)1.0- ((NumericT)1.0-elem_geq0(ul))*elem_geq0(ur) );
    hatf = ustar * ustar / 2 ;
    }

    // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    mesh<NumericT> hatfdiff{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 0) = -1; tKernel(0, 1) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,0}, {0,1}};
    hatfdiff = convolve<NumericT, ConvolutionType::EQUIV>(hatf, tKernel, tROIrc_vec);
    }

    // Time Forward
    return iMesh - hatfdiff * (dt / dx);
}

// TODO  Make it generic for various function 
template <typename NumericT, FluxType fluxT>
inline mesh<NumericT> FluxFunction(const mesh<NumericT> & u_neg, const mesh<NumericT> & u_pos, NumericT dt=0, NumericT dx=0)
{
    if constexpr (fluxT == FluxType::LaxFriedrichs) {
        /* LaxFriedrichs Flux function 
        * \hat{f}_{j+1/2} = \hat{f}(u_{j}, u_{j+1}) 
        * = 1/2 [f(u_j) + f(u_{j+1}) - \alpha (u_{j+1}-u_j)], where \alpha = \max_u |f'(u)|, in the Burger's equation, \alpha = \max(u_j,u_{j+1})
        */
        return (u_neg*u_neg + u_pos*u_pos)/(NumericT)4.0 - elem_max(elem_abs(u_neg), elem_abs(u_pos)) * (u_pos-u_neg);
    } else if constexpr (fluxT == FluxType::LaxWendroff) {
        /* LaxWendroff Flux function 
        * \hat{f}_{j+1/2} = \hat{f}(u_{j}, u_{j+1}) 
        * = 1/2 [f(u_j) + f(u_{j+1})] - dt/2dx f'((u_j+u_{j+1})/2)[f(u_{j+1})-f(u_j)] 
        * Here, due to FVM, u_{j} and u_{j+1} are taken by average value \bar{u}_j and \bar{u}_{j+1}.
        */
        return ((u_neg^2) + (u_pos^2))/4.0 - (dt/dx/2 * 0.5 * 0.5) * (u_neg + u_pos) * ((u_pos^2) - (u_neg^2));
    }
}

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
    // Calculate f(u^*)
    mesh<NumericT> hatf{iMesh.get_size_num()};
    {
    // Finite Volumn 3 order estimation of the cell edge value.
    viennacl::matrix<NumericT> tKernel1{1, 5}; tKernel1(0,1)=-1.0/6; tKernel1(0,2)=5.0/6; tKernel1(0,3) = 1.0/3;
    viennacl::matrix<NumericT> tKernel2{1, 5}; tKernel2(0,2)= 1.0/3; tKernel2(0,3)=5.0/6; tKernel2(0,4) =-1.0/6;
    std::vector<cord2<GridIntT>> tROIrc_vec1{{0,1},{0,2},{0,3}}, tROIrc_vec2{{0,2},{0,3},{0,4}};
    // u_{j+1/2}^- = -1/6 * \bar{u}_{j-1} + 5/6 * \bar{u}_j     + 1/3 * \bar{u}_{j+1}
    mesh<NumericT> u_neg = convolve(iMesh, tKernel1, tROIrc_vec1);
    // u_{j+1/2}^+ =  1/6 * \bar{u}_{j}   + 5/6 * \bar{u}_{j+1} - 1/6 * \bar{u}_{j+2}
    mesh<NumericT> u_pos = convolve(iMesh, tKernel2, tROIrc_vec2);
    
    hatf = FluxFunction<NumericT, FluxType::LaxWendroff>(u_neg, u_pos);
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



enum EdgeFaceOrien : bool { POS=true, NEG=false };

template <typename NumericT, EdgeFaceOrien pos>
/** @brief Finite Volumn 3 order estimation of the cell edge value.
 * 
 * @param  {const mesh<NumericT> &} iMesh : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> FV_EdgeEstimate(const mesh<NumericT> & iMesh)
{
    viennacl::matrix<NumericT> tKernel1{1, 5}; tKernel1(0,1)=-1.0/6; tKernel1(0,2)=5.0/6; tKernel1(0,3) = 1.0/3;
    viennacl::matrix<NumericT> tKernel2{1, 5}; tKernel2(0,2)= 1.0/3; tKernel2(0,3)=5.0/6; tKernel2(0,4) =-1.0/6;
    std::vector<cord2<GridIntT>> tROIrc_vec1{{0,1},{0,2},{0,3}}, tROIrc_vec2{{0,2},{0,3},{0,4}};
    if constexpr (pos == EdgeFaceOrien::NEG) {
        // u_{j+1/2}^- = -1/6 * \bar{u}_{j-1} + 5/6 * \bar{u}_j     + 1/3 * \bar{u}_{j+1}
        return convolve(iMesh, tKernel1, tROIrc_vec1);
    } else {
        // u_{j+1/2}^+ =  1/6 * \bar{u}_{j}   + 5/6 * \bar{u}_{j+1} - 1/6 * \bar{u}_{j+2}
        return convolve(iMesh, tKernel2, tROIrc_vec2);
    }
}

template <typename NumericT, EdgeFaceOrien pos>
/** @brief Finite Volumn 5 order estimation of the cell edge value.
 * 
 * @param  {const mesh<NumericT> &} iMesh : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> WENO_g(const mesh<NumericT> & iMesh)
{

    viennacl::matrix<NumericT> tKernel0{1, 5}; 
    viennacl::matrix<NumericT> tKernel1{1, 5}; 
    viennacl::matrix<NumericT> tKernel2{1, 5}; 
    if constexpr (pos==EdgeFaceOrien::POS) {
        tKernel0(0,4)= 1.0/3; tKernel0(0,3)=-7.0/6; tKernel0(0,2) = 11.0/6;
        tKernel1(0,3)=-1.0/6; tKernel1(0,2)= 5.0/6; tKernel1(0,1) = 1.0/3;
        tKernel2(0,2)= 1.0/3; tKernel2(0,1)= 5.0/6; tKernel2(0,0) =-1.0/6;
    } else {
        tKernel0(0,0)= 1.0/3; tKernel0(0,1)=-7.0/6; tKernel0(0,2) = 11.0/6;
        tKernel1(0,1)=-1.0/6; tKernel1(0,2)= 5.0/6; tKernel1(0,3) = 1.0/3;
        tKernel2(0,2)= 1.0/3; tKernel2(0,3)= 5.0/6; tKernel2(0,4) =-1.0/6;
    }

    static std::vector<cord2<GridIntT>> 
        tROIrc_vec0 = pos ? std::vector<cord2<GridIntT>>{{0,2},{0,3},{0,4}} : std::vector<cord2<GridIntT>>{{0,0},{0,1},{0,2}}, 
        tROIrc_vec1 = std::vector<cord2<GridIntT>>{{0,1},{0,2},{0,3}}, 
        tROIrc_vec2 = pos ? std::vector<cord2<GridIntT>>{{0,0},{0,1},{0,2}} : std::vector<cord2<GridIntT>>{{0,2},{0,3},{0,4}};

    auto u_bracket_0 = convolve(iMesh, tKernel0, tROIrc_vec0); 
    auto u_bracket_1 = convolve(iMesh, tKernel1, tROIrc_vec1);
    auto u_bracket_2 = convolve(iMesh, tKernel2, tROIrc_vec2);

    viennacl::matrix<NumericT> betaKer01{1, 5}; 
    viennacl::matrix<NumericT> betaKer11{1, 5}; 
    viennacl::matrix<NumericT> betaKer21{1, 5}; 
    viennacl::matrix<NumericT> betaKer02{1, 5}; 
    viennacl::matrix<NumericT> betaKer12{1, 5}; 
    viennacl::matrix<NumericT> betaKer22{1, 5}; 
    if constexpr (pos==EdgeFaceOrien::POS) {
        betaKer01(0,4)= 1.0; betaKer01(0,3)=-2.0; betaKer01(0,2) = 1.0;
        betaKer11(0,3)= 1.0; betaKer11(0,2)=-2.0; betaKer11(0,1) = 1.0;
        betaKer21(0,2)= 1.0; betaKer21(0,1)=-2.0; betaKer21(0,0) = 1.0;
        betaKer02(0,4)= 1.0; betaKer02(0,3)=-4.0; betaKer02(0,2) = 3.0;
        betaKer12(0,3)=1.0;/*betaKer12(0,2)=-2.0;*/betaKer12(0,1)=-1.0;
        betaKer22(0,2)= 3.0; betaKer22(0,1)=-4.0; betaKer22(0,0) = 1.0;
    } else {
        betaKer01(0,0)= 1.0; betaKer01(0,1)=-2.0; betaKer01(0,2) = 1.0;
        betaKer11(0,1)= 1.0; betaKer11(0,2)=-2.0; betaKer11(0,3) = 1.0;
        betaKer21(0,2)= 1.0; betaKer21(0,3)=-2.0; betaKer21(0,4) = 1.0;
        betaKer02(0,0)= 1.0; betaKer02(0,1)=-4.0; betaKer02(0,2) = 3.0;
        betaKer12(0,1)=1.0;/*betaKer12(0,2)=-2.0;*/betaKer12(0,3)=-1.0;
        betaKer22(0,2)= 3.0; betaKer22(0,3)=-4.0; betaKer22(0,4) = 1.0;
    }
    static std::vector<cord2<GridIntT>> tROIrc_vec3{{0,1},{0,3}};
    auto beta0_1 = convolve(iMesh, betaKer01, tROIrc_vec0); auto beta0_2 = convolve(iMesh, betaKer02, tROIrc_vec0);
    auto beta1_1 = convolve(iMesh, betaKer11, tROIrc_vec1); auto beta1_2 = convolve(iMesh, betaKer12, tROIrc_vec3);
    auto beta2_1 = convolve(iMesh, betaKer21, tROIrc_vec2); auto beta2_2 = convolve(iMesh, betaKer22, tROIrc_vec2);
    auto beta0 = (NumericT)(13.0/12) * (beta0_1^2) + (NumericT)(1.0/4) * (beta0_2^2);
    auto beta1 = (NumericT)(13.0/12) * (beta1_1^2) + (NumericT)(1.0/4) * (beta1_2^2);
    auto beta2 = (NumericT)(13.0/12) * (beta2_1^2) + (NumericT)(1.0/4) * (beta2_2^2);
    static const NumericT 
        gamma_0 = pos ? 3.0/10 : 1.0/10, 
        gamma_1 = 3.0/5, 
        gamma_2 = pos ? 1.0/10 : 3.0/10,
        epsilon = 1.0e-6;
    auto tilde_w_0 =  gamma_0 / ((beta0+epsilon)^2); 
    auto tilde_w_1 =  gamma_1 / ((beta1+epsilon)^2); 
    auto tilde_w_2 =  gamma_2 / ((beta2+epsilon)^2); 
    {
        auto sum_tilde_w= tilde_w_0 + tilde_w_1 + tilde_w_2;
        tilde_w_0 /= sum_tilde_w; // Now \tilde{w}_0 => w_0 
        tilde_w_1 /= sum_tilde_w;
        tilde_w_2 /= sum_tilde_w;
    }
    if constexpr (pos){ // return u_pos
        viennacl::matrix<NumericT> final_shift_for_pos{1,3}; final_shift_for_pos(0,2)=1.0;
        std::vector<cord2<GridIntT>> final_shift_ROIrc_vec{{0,2}};
        // return tilde_w_0 * u_bracket_0 + tilde_w_1 * u_bracket_1 + tilde_w_2 * u_bracket_2;
        return convolve(
               tilde_w_0 * u_bracket_0 + tilde_w_1 * u_bracket_1 + tilde_w_2 * u_bracket_2, 
               final_shift_for_pos, final_shift_ROIrc_vec);
    } else { // return u_neg
        return tilde_w_0 * u_bracket_0 + tilde_w_1 * u_bracket_1 + tilde_w_2 * u_bracket_2;
    }
}

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
    
    // Minimod function as TVD limiter
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
    u_neg = iMesh          + minimod(u_neg - iMesh,       u_jp1_m_u_j, u_j_m_u_jm1);
    u_pos = convolve(iMesh - minimod(iMesh - u_jm1_pos,   u_jp1_m_u_j, u_j_m_u_jm1), tKernel6, tROIrc_vec6);

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
