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

#include "viennapde/core/cord.hpp"
#include "viennapde/core/mesh.hpp"

namespace viennapde
{

enum ConvolutionType
{
    INNER,
    OUTER,
    EQUIV
};


template < typename NumericT, viennapde::ConvolutionType convolT >
cord2<size_t> ConvolOMatSize(
    const size_t MatSize1, const size_t MatSize2,
    const size_t KerSize1, const size_t KerSize2)
{
    assert( MatSize1>=1 && MatSize2>=1 && KerSize1>=1 && KerSize2>=1 );
    assert( KerSize1 % 2 == 1 && KerSize2 % 2 == 1 );
    const size_t iKernelHalf1 = (KerSize1-1)/2;
    const size_t iKernelHalf2 = (KerSize2-1)/2;
    size_t oMatSize1, oMatSize2;
    if constexpr (convolT==ConvolutionType::EQUIV)
    {
        oMatSize1 = MatSize1;
        oMatSize2 = MatSize2;
    }
    else if constexpr (convolT==ConvolutionType::INNER)
    {
        oMatSize1 = MatSize1 - 2 * iKernelHalf1;
        oMatSize2 = MatSize2 - 2 * iKernelHalf2;
        assert(oMatSize1>=1 && oMatSize2>=1);
    }
    else if constexpr (convolT==ConvolutionType::OUTER)
    {
        oMatSize1 = MatSize1 + 2 * iKernelHalf1;
        oMatSize2 = MatSize2 + 2 * iKernelHalf2;
    }
    return cord2<size_t>(oMatSize1, oMatSize2);
    
}

template < typename NumericT, viennapde::ConvolutionType convolT >
cord2<size_t> ConvolOMatSize(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel)
{
    return ConvolOMatSize<NumericT, convolT>(
        iMatrix.size1(), iMatrix.size2(), 
        iKernel.size1(), iKernel.size2());
}

template < typename NumericT, viennapde::ConvolutionType convolT >
cord3<size_t> ConvolOMeshSize(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennapde::Varmesh<NumericT> & iKernel)
{
    const cord2<size_t> xySize= ConvolOMatSize<NumericT, convolT>(*iVarmesh[0], *iKernel[0]);
    assert( iKernel.size() % 2 == 1 );
    const size_t MeshSize3 = iKernel.size(),
                   iKernelHalf3 = (MeshSize3-1)/2;
          size_t oMeshSize3;
    if constexpr (convolT==ConvolutionType::EQUIV)
    {
        oMeshSize3 = MeshSize3;
    }
    else if constexpr (convolT==ConvolutionType::INNER)
    {
        oMeshSize3 = MeshSize3 - 2 * iKernelHalf3;
        assert(oMeshSize3>=1);
    }
    else if constexpr (convolT==ConvolutionType::OUTER)
    {
        oMeshSize3 = MeshSize3 + 2 * iKernelHalf3;
    }
    return cord3<size_t>{xySize.x, xySize.y, (size_t)oMeshSize3}; 
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
    std::vector<cord2<GridIntT>>& ROIrc_vec, ClrOut clrOut = ClrOut::YES) // TODO ROIrc_vec, I am worrying that for the full kernel, this additional option may be speed bottleneck. 
{
    assert( &iMatrix != &oMatrix ); // TODO: This function not yet designed for self convolving.
    // STUB 01 Check the matrix size 
    assert( (iKernel.size1() % 2 == 1) && (iKernel.size2() % 2 == 1) );
    const size_t    iKernelHalf1 = (iKernel.size1()-1)/2,
                    iKernelHalf2 = (iKernel.size2()-1)/2;
    if (clrOut==ClrOut::YES) {
        if constexpr (convolT==EQUIV) {
            oMatrix.resize(iMatrix.size1(), iMatrix.size2(), false);
        } else if constexpr (convolT==INNER) {
            oMatrix.resize(iMatrix.size1()+2*iKernelHalf1, iMatrix.size2()+2*iKernelHalf2, false);
        } else if constexpr (convolT==OUTER) {
            oMatrix.resize(iMatrix.size1()-2*iKernelHalf1, iMatrix.size2()-2*iKernelHalf2, false);
        }
    } else {
        if constexpr (convolT==EQUIV) {
            assert(iMatrix.size1()==oMatrix.size1() && iMatrix.size2()==oMatrix.size2());
        } else if constexpr (convolT==INNER) {
            assert(iMatrix.size1()==oMatrix.size1()+2*iKernelHalf1 && iMatrix.size2()==oMatrix.size2()+2*iKernelHalf2);
        } else if constexpr (convolT==OUTER) {
            assert(iMatrix.size1()==oMatrix.size1()-2*iKernelHalf1 && iMatrix.size2()==oMatrix.size2()-2*iKernelHalf2);
        }
    }

    
    
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.
    if (clrOut) {oMatrix.clear();}

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
viennacl::matrix<NumericT> convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<cord2<GridIntT>>& ROIrc_vec, ClrOut clrOut = ClrOut::YES) //TODO: I guess the TODO marker can be removed. 
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(iMatrix, iKernel);    
    viennacl::matrix<NumericT> oMatrix{(size_t)oMatSize.x, (size_t)oMatSize.y};
    viennapde::convolve<NumericT, convolT>(iMatrix, iKernel, oMatrix, ROIrc_vec, clrOut);
    return oMatrix;
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
void convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel,
    viennacl::matrix<NumericT> & oMatrix, ClrOut clrOut = ClrOut::YES) 
{
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve<NumericT, convolT>(iMatrix, iKernel, oMatrix, ROIrc_vec, clrOut);
} //function void viennapde::convolve

template <  typename NumericT, 
            viennapde::ConvolutionType convolT = EQUIV>
viennacl::matrix<NumericT> convolve(
    const viennacl::matrix<NumericT> & iMatrix,
    const viennacl::matrix<NumericT> & iKernel) 
{
    const cord2<size_t> oMatSize = ConvolOMatSize<NumericT, convolT>(iMatrix, iKernel);    
    viennacl::matrix<NumericT> oMatrix{oMatSize.x, oMatSize.y};
    std::vector<cord2<GridIntT>> ROIrc_vec{}; 
    for (size_t i = 0; i < iKernel.size1(); i++)
    for (size_t j = 0; j < iKernel.size2(); j++)
        ROIrc_vec.push_back(cord2<GridIntT>(i, j));
    viennapde::convolve<NumericT, convolT>(iMatrix, iKernel, oMatrix, ROIrc_vec);
    return oMatrix;
} //function void viennapde::convolve


} //namespace viennapde


