#pragma once
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

/** @file viennapde/core/convol_mesh.hpp
    @brief Convolve operation on Varmesh class. 
    All functions in this file are overloaded of viennapde::convolve working on the Varmesh class. 
    TODO Not yet and verified.
*/
#include <vector>

#include "viennacl/matrix.hpp"
#include "./mesh.hpp"
#include "./convol_mat.hpp"

namespace viennapde
{

/*===== SECTION Convolve Mesh --Mat-> Mesh ==================================================== */
// TODO totally change the comment & plut the output initialization & temp Varmesh removal
// SECTION 03_002a Varmesh Convolution
/** @brief Convolve the Varmesh data by the 2D matrix kernel, which would be the base of Varmesh shift and filter
 * 
 * @param  {viennacl::matrix<NumericT>} iKernel    : 
 * @param  {std::vector<std::pair<size_t} undefined : 
 * @param  {size_t>>} ROIxy_vec                     : 
 * 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh,
    std::vector<cord2<GridIntT>> & ROIrc_vec, ClrOut clrOut = ClrOut::YES)
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(*iVarmesh[0], iKernel);
    oVarmesh.resize_ptr(iVarmesh.get_layer_num(), oMatSize.x, oMatSize.y); // I suggest to use a completely new mesh object.
    if (clrOut == ClrOut::YES)
    {
        for (size_t layer_i=0; layer_i< iVarmesh.get_layer_num(); layer_i++) 
            *(oVarmesh[layer_i]) = viennapde::convolve<NumericT, convolT>(*(iVarmesh[layer_i]), iKernel, ROIrc_vec, clrOut);
    } else //(clrOut == ClrOut::NO)
    {
        for (size_t layer_i=0; layer_i< iVarmesh.get_layer_num(); layer_i++) 
            *(oVarmesh[layer_i]) += viennapde::convolve<NumericT, convolT>(*(iVarmesh[layer_i]), iKernel, ROIrc_vec, clrOut);
    }
    
} //function void viennapde::convolve



// SECTION 01_002b Mesh Convolution
/** @brief the same as the convolve funciton 01_002a, just the default ROI becomes the whole matrix.
 * @param  {viennacl::matrix<NumericT>} iKernel : 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::Varmesh<NumericT> convolve(
    viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<cord2<GridIntT>> & ROIrc_vec)
{
    viennapde::Varmesh<NumericT> tVarmesh(iVarmesh.get_size_num());
    viennapde::convolve<NumericT, convolT>(iVarmesh, iKernel, tVarmesh, ROIrc_vec);
    return tVarmesh;
} //function viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh, ClrOut clrOut = ClrOut::YES) 
{
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iVarmesh, iKernel, oVarmesh, ROIrc_vec, clrOut);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::Varmesh<NumericT> convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel) 
{
    const std::pair <size_t, size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(
        iVarmesh.get_row_num(), iVarmesh.get_column_num(), iKernel.size1(), iKernel.size2());    
    viennacl::matrix<NumericT> oVarmesh{oMatSize.first, oMatSize.second, iVarmesh.get_layer_num()};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iVarmesh, iKernel, oVarmesh, ROIrc_vec, ClrOut::YES);
    return oVarmesh;
} //function void viennapde::convolve



/*===== SECTION Convolve Mesh --Mesh-> Mesh ==================================================== */

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennapde::Varmesh<NumericT> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh,
    std::vector<cord2<GridIntT>> & ROIrc_vec, ClrOut clrOut = ClrOut::YES)
{
    cord3<size_t> oMeshSize = ConvolOMeshSize(iVarmesh, iKernel);
    oVarmesh.resize_ptr(oMeshSize.z, oMeshSize.x, oMeshSize.y);
    if (clrOut == ClrOut::YES) oVarmesh.clear();
    if constexpr (convolT==ConvolutionType::OUTER)
    {
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
            viennapde::convolve<NumericT, convolT>(
                iVarmesh, 
                iKernel[layer_i+(iKernel.get_layer_num()-1)/2], 
                viennapde::Varmesh<NumericT>{
                    oVarmesh.begin() + (iKernel.get_layer_num()-1)/2 - layer_i, 
                    oVarmesh.begin() + (iKernel.get_layer_num()-1)/2 - layer_i + iVarmesh.get_layer_num()}, 
                ROIrc_vec, ClrOut::NO);
    } else if constexpr (convolT==ConvolutionType::EQUIV)
    {   
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
            viennapde::convolve<NumericT, convolT>(
                viennapde::Varmesh<NumericT>{
                    iVarmesh.begin() + std::max(layer_i,(ssize_t)0), 
                    iVarmesh.begin() + std::min(layer_i,(ssize_t)0) + iVarmesh.get_layer_num()}, 
                iKernel[layer_i+(iKernel.get_layer_num()-1)/2], 
                viennapde::Varmesh<NumericT>{
                    oVarmesh.begin() - std::min(layer_i,(ssize_t)0), 
                    oVarmesh.begin() - std::max(layer_i,(ssize_t)0) + iVarmesh.get_layer_num()}, 
                ROIrc_vec, ClrOut::NO);
    } else if constexpr (convolT==ConvolutionType::INNER)
    {   
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
            viennapde::convolve<NumericT, convolT>(
                viennapde::Varmesh<NumericT>{
                    iVarmesh.begin() + layer_i + (iKernel.get_layer_num()-1)/2, 
                    iVarmesh.begin() + layer_i + (iKernel.get_layer_num()-1)/2 + oVarmesh.get_layer_num()}, 
                iKernel[layer_i+(iKernel.get_layer_num()-1)/2], 
                oVarmesh, 
                ROIrc_vec, ClrOut::NO);
    }
    
     
} //function void viennapde::convolve


template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::Varmesh<NumericT> convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennapde::Varmesh<NumericT> & iKernel,
    std::vector<cord2<GridIntT>> & ROIrc_vec)
{
    viennapde::Varmesh<NumericT> tVarmesh{};
    viennapde::convolve<NumericT, convolT>(iVarmesh, iKernel, tVarmesh, ROIrc_vec);
    return tVarmesh;
}


// TODO The Mesh convolved by mesh operation is contains kernel elements as 3D cord. Change the cord from 2D to 3D.
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennapde::Varmesh<NumericT> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh, ClrOut clrOut = ClrOut::YES) 
{
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iVarmesh, iKernel, oVarmesh, ROIrc_vec, clrOut);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::Varmesh<NumericT> convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennapde::Varmesh<NumericT> & iKernel) 
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(
        iVarmesh.get_row_num(), iVarmesh.get_column_num(), iKernel.size1(), iKernel.size2());    
    viennacl::matrix<NumericT> oVarmesh{oMatSize.x, oMatSize.y, iVarmesh.get_layer_num()};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iVarmesh, iKernel, oVarmesh, ROIrc_vec, ClrOut::YES);
    return oVarmesh;
} //function void viennapde::convolve

} //namespace viennapde


