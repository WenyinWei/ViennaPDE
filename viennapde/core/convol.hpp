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






// SECTION 03_002a Vermesh Convolution, utilizing 01_002a
/** @brief Convolve the Varmesh data by the 2D matrix kernel, which would be the base of Varmesh shift and filter
 * @param  {viennacl::matrix<NumericT>} iKernel    : 
 * @param  {std::vector<std::pair<size_t} undefined : 
 * @param  {size_t>>} ROIxy_vec                     : 
 * 
 * @example
 * vcl_MatrixT kernel(5, 5);
 * kernel.clear();
 * kernel(2, 1) = -1.0/2; kernel(2, 2) = 0; kernel(2, 3) = 1.0/2;
 * std::vector<std::pair<size_t, size_t>> vec;
 * vec.push_back(std::make_pair<int,int>(2, 3));
 * vcl_varmesh.convolve(kernel, vec);
 */
template <  typename NumericT, 
            viennapde::ConvolutionType ConvolType = EQUIV>
void convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    viennacl::matrix<NumericT> & oMatrix,
    std::vector<cord2<GridIntT>> ROIrc_vec = std::vector<cord2<GridIntT>>() ) 
{
    // STUB 01 
    assert( (iKernel.size1() % 2 == 1) && (iKernel.size2() % 2 == 1) );
    const size_t iKernelHalf1 = (iKernel.size1()-1)/2;
    const size_t iKernelHalf2 = (iKernel.size2()-1)/2;
    if constexpr (ConvolType==EQUIV)
    {
        assert(iMatrix.size1()==oMatrix.size1());
        assert(iMatrix.size2()==oMatrix.size2());
    }

    // NOTE Argument ROIrc_vec default empty case, all entries are filled in it here.
    if (ROIrc_vec.empty())
    {
        for (size_t i = 0; i < iKernel.size1(); i++)
        for (size_t j = 0; j < iKernel.size2(); j++)
            ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    }
    
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.

    oMatrix.clear(); // REVIEW This may be the efficiency bottleneck which is safe but slow 
    for(auto & iter: ROIrc_vec)
    {
        size_t  l_row = iter.x, 
                l_column = iter.y;
        int bias1 = l_row-iKernelHalf1; 
        int bias2 = l_column-iKernelHalf2;
        // REVIEW The following command is helpful to avoid unexpected too big kernel, however, I am not satisfactory with the if statement which may be a bottleneck in this loop.
        if ((std::abs(bias1)>=iKernelHalf1) || (std::abs(bias2)>=iKernelHalf2)) continue;
        if constexpr (ConvolType == ConvolutionType::EQUIV)
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
        else // TODO: Not yet implemented for other convolution type
        {
            std::cerr << "Not yet implemented for other convolution type than EQUIV." << std::endl;
            viennacl::range submat_row_range(l_row, l_row+iMatrix.size1()), 
                            submat_column_range(l_column, l_column+iMatrix.size2());
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(iMatrix, submat_row_range, submat_column_range); 
            oMatrix += iKernel(l_row, l_column) * submatrix_to_add;
        }
        
    }    
} //function void viennapde::convolve




} //namespace viennapde


