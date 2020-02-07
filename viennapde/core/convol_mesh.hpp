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
    // TODO take care memory newing
    if (oVarmesh.get_size_num() == iVarmesh.get_size_num())
    {
        /* code */
    }
    
    oVarmesh.resize_ptr(iVarmesh.get_layer_num()); 
    // FIXME I am worried here that the shared object is not created properly.
    for (size_t layer_i=0; layer_i< iVarmesh.get_layer_num(); layer_i++) 
        *(oVarmesh[layer_i])
        = viennapde::convolve<NumericT, convolT>(*(iVarmesh[layer_i]),iKernel, ROIrc_vec, clrOut);
} //function void viennapde::convolve



// SECTION 01_002b Mesh Convolution
/** @brief the same as the convolve funciton 01_002a, just the default ROI becomes the whole matrix.
 * @param  {viennacl::matrix<NumericT>} iKernel : 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennapde::Varmesh<NumericT> & convolve(
    viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<cord2<GridIntT>> & ROIrc_vec, ClrOut clrOut = ClrOut::YES)
{
    viennapde::Varmesh<NumericT> * tVarmesh(iVarmesh.get_layer_num(),iVarmesh.get_row_num(), iVarmesh.get_column_num());
    viennapde::convolve<NumericT, convolT>(iVarmesh, iKernel, tVarmesh, ROIrc_vec, clrOut);
    return *tVarmesh;
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
viennapde::Varmesh<NumericT> & convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel) 
{
    const std::pair <size_t, size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(
        iVarmesh.get_row_num(), iVarmesh.get_column_num(), iKernel.size1(), iKernel.size2());    
    viennacl::matrix<NumericT> *oVarmesh = new viennapde::Varmesh<NumericT> {iVarmesh.get_layer_num(), oMatSize.first, oMatSize.second};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iVarmesh, iKernel, *oVarmesh, ROIrc_vec, ClrOut::YES);
    return *oVarmesh;
} //function void viennapde::convolve



// TODO The following functions are for 3D mesh convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const std::deque<viennacl::matrix<NumericT>> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh,
    std::vector<cord2<GridIntT>> & ROIrc_vec, ClrOut clrOut = ClrOut::YES)
{
    assert( iKernel.size() % 2 == 1 );
    oVarmesh.data_->resize(iVarmesh.get_layer_num()); 
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.
    for (size_t layer_i=0; layer_i< iVarmesh.get_layer_num(); layer_i++) 
        oVarmesh.data_->at(layer_i) 
        = viennapde::convolve<NumericT, convolT>(iVarmesh.data_->at(layer_i),iKernel, ROIrc_vec, clrOut);
} //function void viennapde::convolve

} //namespace viennapde


