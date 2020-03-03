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

/** @file viennapde/scheme/pde/edge_value_estimate.hpp
    @brief Methods based on FV or FD method to estimate the conservative value at the cell edge (including both -/+ faces).
*/

#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mesh.hpp"


// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

enum EdgeFaceOrien : bool { POS=true, NEG=false };

template <typename NumericT, EdgeFaceOrien edge_face_orien>
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
    if constexpr (edge_face_orien == EdgeFaceOrien::NEG) {
        // u_{j+1/2}^- = -1/6 * \bar{u}_{j-1} + 5/6 * \bar{u}_j     + 1/3 * \bar{u}_{j+1}
        return convolve(iMesh, tKernel1, tROIrc_vec1);
    } else {
        // u_{j+1/2}^+ =  1/6 * \bar{u}_{j}   + 5/6 * \bar{u}_{j+1} - 1/6 * \bar{u}_{j+2}
        return convolve(iMesh, tKernel2, tROIrc_vec2);
    }
}

template <typename NumericT, EdgeFaceOrien edge_face_orien>
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
    if constexpr (edge_face_orien==EdgeFaceOrien::POS) {
        tKernel0(0,4)= 1.0/3; tKernel0(0,3)=-7.0/6; tKernel0(0,2) = 11.0/6;
        tKernel1(0,3)=-1.0/6; tKernel1(0,2)= 5.0/6; tKernel1(0,1) = 1.0/3;
        tKernel2(0,2)= 1.0/3; tKernel2(0,1)= 5.0/6; tKernel2(0,0) =-1.0/6;
    } else {
        tKernel0(0,0)= 1.0/3; tKernel0(0,1)=-7.0/6; tKernel0(0,2) = 11.0/6;
        tKernel1(0,1)=-1.0/6; tKernel1(0,2)= 5.0/6; tKernel1(0,3) = 1.0/3;
        tKernel2(0,2)= 1.0/3; tKernel2(0,3)= 5.0/6; tKernel2(0,4) =-1.0/6;
    }

    static std::vector<cord2<GridIntT>> 
        tROIrc_vec0 = edge_face_orien ? std::vector<cord2<GridIntT>>{{0,2},{0,3},{0,4}} : std::vector<cord2<GridIntT>>{{0,0},{0,1},{0,2}}, 
        tROIrc_vec1 = std::vector<cord2<GridIntT>>{{0,1},{0,2},{0,3}}, 
        tROIrc_vec2 = edge_face_orien ? std::vector<cord2<GridIntT>>{{0,0},{0,1},{0,2}} : std::vector<cord2<GridIntT>>{{0,2},{0,3},{0,4}};

    auto u_bracket_0 = convolve(iMesh, tKernel0, tROIrc_vec0); 
    auto u_bracket_1 = convolve(iMesh, tKernel1, tROIrc_vec1);
    auto u_bracket_2 = convolve(iMesh, tKernel2, tROIrc_vec2);

    viennacl::matrix<NumericT> betaKer01{1, 5}; 
    viennacl::matrix<NumericT> betaKer11{1, 5}; 
    viennacl::matrix<NumericT> betaKer21{1, 5}; 
    viennacl::matrix<NumericT> betaKer02{1, 5}; 
    viennacl::matrix<NumericT> betaKer12{1, 5}; 
    viennacl::matrix<NumericT> betaKer22{1, 5}; 
    if constexpr (edge_face_orien==EdgeFaceOrien::POS) {
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
    auto beta0 = (13.0/12) * (beta0_1^2) + (1.0/4) * (beta0_2^2);
    auto beta1 = (13.0/12) * (beta1_1^2) + (1.0/4) * (beta1_2^2);
    auto beta2 = (13.0/12) * (beta2_1^2) + (1.0/4) * (beta2_2^2);
    static const NumericT 
        gamma_0 = edge_face_orien ? 3.0/10 : 1.0/10, 
        gamma_1 = 3.0/5, 
        gamma_2 = edge_face_orien ? 1.0/10 : 3.0/10,
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
    if constexpr (edge_face_orien){ // return u_pos
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

} //namespace viennapde::scheme
} //namespace viennapde
