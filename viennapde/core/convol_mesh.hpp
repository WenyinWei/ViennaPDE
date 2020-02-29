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
    @brief Convolve operation on mesh class. 
    All functions in this file are overloaded of viennapde::convolve working on the mesh class. 
    TODO Not yet and verified.
*/
#include <vector>

#include "viennacl/matrix.hpp"
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mat.hpp"

namespace viennapde
{

/*===== SECTION Convolve Mat --on-- Mat-> Mesh ==================================================== */
// TODO totally change the comment & plut the output initialization & temp mesh removal
// SECTION 03_002a mesh Convolution
/** @brief Convolve the mesh data by the 2D matrix kernel, which would be the base of mesh shift and filter
 * 
 * @param  {viennacl::matrix<NumericT>} iKernel    : 
 * @param  {std::vector<std::pair<size_t} undefined : 
 * @param  {size_t>>} ROIxy_vec                     : 
 * 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennacl::matrix<NumericT> & iKernel,
    viennapde::mesh<NumericT> & oMesh,
    std::vector<cord2<GridIntT>> & ROIrc_vec, ClrOut clrOut = ClrOut::YES)
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(*(iMesh[0]), iKernel);
    oMesh.resize_ptr(iMesh.get_layer_num(), oMatSize.x, oMatSize.y); // I suggest to use a completely new mesh object.
    if (clrOut == ClrOut::YES)
    {
        for (size_t layer_i=0; layer_i< iMesh.get_layer_num(); layer_i++) 
            *(oMesh[layer_i]) = viennapde::convolve<NumericT, convolT>(*(iMesh[layer_i]), iKernel, ROIrc_vec, clrOut);
    } else //(clrOut == ClrOut::NO)
    {
        for (size_t layer_i=0; layer_i< iMesh.get_layer_num(); layer_i++) 
            *(oMesh[layer_i]) += viennapde::convolve<NumericT, convolT>(*(iMesh[layer_i]), iKernel, ROIrc_vec, clrOut);
    }
    
} //function void viennapde::convolve



// SECTION 01_002b Mesh Convolution
/** @brief the same as the convolve funciton 01_002a, just the default ROI becomes the whole matrix.
 * @param  {viennacl::matrix<NumericT>} iKernel : 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::mesh<NumericT> convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<cord2<GridIntT>> & ROIrc_vec)
{
    viennapde::mesh<NumericT> tMesh(iMesh.get_size_num());
    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, tMesh, ROIrc_vec);
    return tMesh;
} //function viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennacl::matrix<NumericT> & iKernel,
    viennapde::mesh<NumericT> & oMesh, ClrOut clrOut = ClrOut::YES) 
{
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.emplace_back(i, j);
    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, oMesh, ROIrc_vec, clrOut);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::mesh<NumericT> convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennacl::matrix<NumericT> & iKernel) 
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(
        iMesh.get_row_num(), iMesh.get_column_num(), iKernel.size1(), iKernel.size2());    
    viennapde::mesh<NumericT> oMesh{oMatSize.x, oMatSize.y, iMesh.get_layer_num()};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.emplace_back(i, j);
    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, oMesh, ROIrc_vec, ClrOut::YES);
    return oMesh;
} //function void viennapde::convolve



/*===== SECTION Convolve Mesh --on-- Mesh-> Mesh ==================================================== */

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennapde::mesh<NumericT> & iKernel,
    viennapde::mesh<NumericT> & oMesh,
    std::vector<cord3<GridIntT>> & ROIrcl_vec, ClrOut clrOut = ClrOut::YES)
{
    cord3<size_t> oMeshSize = ConvolOMeshSize<NumericT, convolT>(iMesh, iKernel);
    if (clrOut == ClrOut::YES)
    {
        oMesh.resize_ptr(0);
        oMesh.resize_ptr(oMeshSize.z, oMeshSize.x, oMeshSize.y);
    } else // clrOut == ClrOut::NO
    {
        bool sizeMatch = oMesh.get_row_num() == oMeshSize.x 
                        && oMesh.get_column_num()==oMeshSize.y 
                        && oMesh.get_layer_num()==oMeshSize.z;
        if (!sizeMatch) {
            std::cerr << "You can not do the convol operation and accumulate the result on the oMesh because the two mesh size do not match. \n";
        } // If size does not match, then no need to resize the mesh.
    }
    if constexpr (convolT==ConvolutionType::OUTER)
    {
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
        {
            std::vector<cord2<GridIntT>> ROIrc_vec;
            for (size_t i = 0; i < ROIrcl_vec.size(); i++)
                if (ROIrcl_vec[i].z == layer_i+(iKernel.get_layer_num()-1)/2)  
                    ROIrc_vec.emplace_back(ROIrcl_vec[i].x, ROIrcl_vec[i].y);
            viennapde::mesh<NumericT> oVarmesh_part{
                    oMesh.begin() + (iKernel.get_layer_num()-1)/2 - layer_i, 
                    oMesh.begin() + (iKernel.get_layer_num()-1)/2 - layer_i + iMesh.get_layer_num()};
            viennapde::convolve<NumericT, convolT>(
                iMesh, 
                *(iKernel[layer_i+(iKernel.get_layer_num()-1)/2]), 
                oVarmesh_part,
                ROIrc_vec, ClrOut::NO);
        }
    } else if constexpr (convolT==ConvolutionType::EQUIV)
    {   
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
        {
            std::vector<cord2<GridIntT>> ROIrc_vec;
            for (size_t i = 0; i < ROIrcl_vec.size(); i++)
                if (ROIrcl_vec[i].z == layer_i+(iKernel.get_layer_num()-1)/2)  
                    ROIrc_vec.emplace_back(ROIrcl_vec[i].x, ROIrcl_vec[i].y);
            viennapde::mesh<NumericT> iVarmesh_part{
                    iMesh.begin() + std::max(layer_i,(ssize_t)0), 
                    iMesh.begin() + std::min(layer_i,(ssize_t)0) + iMesh.get_layer_num()};
            viennapde::mesh<NumericT> oVarmesh_part{
                    oMesh.begin() - std::min(layer_i,(ssize_t)0), 
                    oMesh.begin() - std::max(layer_i,(ssize_t)0) + iMesh.get_layer_num()};
            viennapde::convolve<NumericT, convolT>(
                iVarmesh_part, 
                *(iKernel[layer_i+(iKernel.get_layer_num()-1)/2]), 
                oVarmesh_part,
                ROIrc_vec, ClrOut::NO);
        }
    } else if constexpr (convolT==ConvolutionType::INNER)
    {   
        for (ssize_t layer_i=-(iKernel.get_layer_num()-1)/2; layer_i< (iKernel.get_layer_num()-1)/2+1; layer_i++) 
        {
            std::vector<cord2<GridIntT>> ROIrc_vec;
            for (size_t i = 0; i < ROIrcl_vec.size(); i++)
                if (ROIrcl_vec[i].z == layer_i+(iKernel.get_layer_num()-1)/2)  
                    ROIrc_vec.emplace_back(ROIrcl_vec[i].x, ROIrcl_vec[i].y);
            viennapde::mesh<NumericT> iVarmesh_part{
                    iMesh.begin() + layer_i + (iKernel.get_layer_num()-1)/2, 
                    iMesh.begin() + layer_i + (iKernel.get_layer_num()-1)/2 + oMesh.get_layer_num()};
            viennapde::convolve<NumericT, convolT>(
                iVarmesh_part, 
                *(iKernel[layer_i+(iKernel.get_layer_num()-1)/2]), 
                oMesh, 
                ROIrc_vec, ClrOut::NO);
        }
    }
    
     
} //function void viennapde::convolve


template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::mesh<NumericT> convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennapde::mesh<NumericT> & iKernel,
    std::vector<cord3<GridIntT>> & ROIrcl_vec)
{
    viennapde::mesh<NumericT> tMesh{};
    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, tMesh, ROIrcl_vec);
    return tMesh;
}


// TODO The Mesh convolved by mesh operation is contains kernel elements as 3D cord. Change the cord from 2D to 3D.
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennapde::mesh<NumericT> & iKernel,
    viennapde::mesh<NumericT> & oMesh, ClrOut clrOut = ClrOut::YES) 
{
    std::vector<cord3<GridIntT>> ROIrcl_vec; 
    for (size_t i = 0; i < iKernel.get_row_num(); i++)
    for (size_t j = 0; j < iKernel.get_column_num(); j++)
    for (size_t k = 0; k < iKernel.get_layer_num(); k++)
        ROIrcl_vec.emplace_back(i, j, k);

    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, oMesh, ROIrcl_vec, clrOut);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::mesh<NumericT> convolve(
    const viennapde::mesh<NumericT> & iMesh,
    const viennapde::mesh<NumericT> & iKernel) 
{
    cord3<size_t> oMeshSize = ConvolOMeshSize<NumericT, convolT>(iMesh, iKernel);
    viennapde::mesh<NumericT> oMesh{oMeshSize.x, oMeshSize.y, oMeshSize.z};
    std::vector<cord3<GridIntT>> ROIrcl_vec{}; 
    for (size_t i = 0; i < iKernel.get_row_num(); i++)
    for (size_t j = 0; j < iKernel.get_column_num(); j++)
    for (size_t k = 0; k < iKernel.get_layer_num(); k++)
        ROIrcl_vec.emplace_back(i, j, k);
    viennapde::convolve<NumericT, convolT>(iMesh, iKernel, oMesh, ROIrcl_vec, ClrOut::YES);
    return oMesh;
} //function void viennapde::convolve

} //namespace viennapde


