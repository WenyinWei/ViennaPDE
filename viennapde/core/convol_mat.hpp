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

/** @file viennapde/core/convol.hpp
    @brief Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
    TODO These functions have not yet been verified.
*/

#include "viennacl/matrix.hpp"



namespace viennapde
{

enum ConvolutionType
{
    INNER,
    OUTER,
    EQUIV
};

template < typename NumericT, viennapde::ConvolutionType convolT >
std::pair<size_t, size_t> ConvolOMatSize(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel)
{
    const size_t iKernelHalf1 = (iKernel.size1()-1)/2;
    const size_t iKernelHalf2 = (iKernel.size2()-1)/2;
    size_t oMatSize1, oMatSize2;
    if constexpr (convolT==ConvolutionType::EQUIV)
    {
        oMatSize1 = iMatrix.size1();
        oMatSize2 = iMatrix.size2();
    }
    else if constexpr (convolT==ConvolutionType::INNER)
    {
        oMatSize1 = iMatrix.size1() - 2 * iKernelHalf1;
        oMatSize2 = iMatrix.size2() - 2 * iKernelHalf2;
    }
    else if constexpr (convolT==ConvolutionType::OUTER)
    {
        oMatSize1 = iMatrix.size1() + 2 * iKernelHalf1;
        oMatSize2 = iMatrix.size2() + 2 * iKernelHalf2;
    }
    return std::pair<size_t, size_t>(oMatSize1, oMatSize2);
    
}

// SECTION 03_002a Vermesh Convolution
/** @brief Convolve the Varmesh data by the 2D matrix kernel, which would be the base of Varmesh shift and filter
 * @param  {viennacl::matrix<NumericT>} iKernel    : 
 * @param  {std::vector<std::pair<size_t} undefined : 
 * @param  {size_t>>} ROIxy_vec                     : 
 */
template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    viennacl::matrix<NumericT> & oMatrix,
    std::vector<cord2<GridIntT>>& ROIrc_vec, bool clearResult = true) // TODO ROIrc_vec, I am worrying that for the full kernel, this additional option may be speed bottleneck. 
{
    // STUB 01 Check the matrix size 
    assert( (iKernel.size1() % 2 == 1) && (iKernel.size2() % 2 == 1) );
    const size_t    iKernelHalf1 = (iKernel.size1()-1)/2,
                    iKernelHalf2 = (iKernel.size2()-1)/2;
    if constexpr (convolT==EQUIV) {
        assert(iMatrix.size1()==oMatrix.size1() && iMatrix.size2()==oMatrix.size2());
    } else if constexpr (convolT==INNER) {
        assert(iMatrix.size1()==oMatrix.size1()+2*iKernelHalf1 && iMatrix.size2()==oMatrix.size2()+2*iKernelHalf2);
    } else if constexpr (convolT==OUTER) {
        assert(iMatrix.size1()==oMatrix.size1()-2*iKernelHalf1 && iMatrix.size2()==oMatrix.size2()-2*iKernelHalf2);
    }
    
    
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.
    if (clearResult) {oMatrix.clear();}

    for(auto & iter: ROIrc_vec)
    {
        size_t  l_row = iter.x, 
                l_column = iter.y;
        int bias1 = l_row-iKernelHalf1; 
        int bias2 = l_column-iKernelHalf2;
        // REVIEW The following command is helpful to avoid unexpected too big kernel, however, I am not satisfactory with the if statement which may be a bottleneck in this loop.
        if ((std::abs(bias1) > iKernelHalf1) || (std::abs(bias2) > iKernelHalf2)) continue;
        if constexpr (convolT == ConvolutionType::EQUIV)
        {
            viennacl::range submat_row_range_from(std::max(bias1, 0), 
                                                    iMatrix.size1() + std::min(bias1, 0)),
                            submat_column_range_from(std::max(bias2, 0), 
                                                    iMatrix.size2() + std::min(bias2, 0)),
                            submat_row_range_to(std::max(-bias1, 0), 
                                                    iMatrix.size1() + std::min(-bias1, 0)),
                            submat_column_range_to(std::max(-bias2, 0), 
                                                    iMatrix.size2() + std::min(-bias2, 0));
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(iMatrix, submat_row_range_from, submat_column_range_from),
                submatrix_tobe_add(oMatrix, submat_row_range_to, submat_column_range_to); 

            submatrix_tobe_add += iKernel(l_row, l_column) * submatrix_to_add;
        }
        else if constexpr (convolT == ConvolutionType::INNER)
        {
            viennacl::range submat_row_range_from(iKernelHalf1+bias1, iKernelHalf1+bias1+oMatrix.size1()),
                            submat_column_range_from(iKernelHalf2+bias2, iKernelHalf2+bias2+oMatrix.size2()),
                            submat_row_range_to(0, oMatrix.size1()),
                            submat_column_range_to(0, oMatrix.size2());
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(iMatrix, submat_row_range_from, submat_column_range_from),
                submatrix_tobe_add(oMatrix, submat_row_range_to, submat_column_range_to); 

            submatrix_tobe_add += iKernel(l_row, l_column) * submatrix_to_add;
        }
        else if constexpr (convolT == ConvolutionType::OUTER)
        {
            viennacl::range submat_row_range_from(0, iMatrix.size1()),
                            submat_column_range_from(0, iMatrix.size2()),
                            submat_row_range_to(iKernelHalf1-bias1, iKernelHalf1+iMatrix.size1()-bias1),
                            submat_column_range_to(iKernelHalf2-bias2, iKernelHalf2+iMatrix.size2()-bias2);
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(iMatrix, submat_row_range_from, submat_column_range_from),
                submatrix_tobe_add(oMatrix, submat_row_range_to, submat_column_range_to); 

            submatrix_tobe_add += iKernel(l_row, l_column) * submatrix_to_add;
        }
        
    }    
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennacl::matrix<NumericT> & convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<cord2<GridIntT>>& ROIrc_vec, bool clearResult = true) //TODO: I guess the TODO marker can be removed. 
{
    const std::pair <size_t, size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(iMatrix, iKernel);    
    viennacl::matrix<NumericT> *oMatrix = new viennacl::matrix<NumericT> {oMatSize.first, oMatSize.second};
    viennapde::convolve<NumericT, convolT>(iMatrix, iKernel, *oMatrix, ROIrc_vec, clearResult);
    return *oMatrix;
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    viennacl::matrix<NumericT> & oMatrix, bool clearResult = true) 
{
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iMatrix, iKernel, oMatrix, ROIrc_vec, clearResult);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennacl::matrix<NumericT> & convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel) 
{
    const std::pair <size_t, size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(iMatrix, iKernel);    
    viennacl::matrix<NumericT> *oMatrix = new viennacl::matrix<NumericT> {oMatSize.first, oMatSize.second};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve(iMatrix, iKernel, *oMatrix, ROIrc_vec);
    return *oMatrix;
} //function void viennapde::convolve


} //namespace viennapde

