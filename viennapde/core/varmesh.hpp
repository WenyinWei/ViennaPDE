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

/** @file viennapde/core/varmesh.hpp
    @brief Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
*/

#include "viennacl/forwards.h"
// #include "viennacl/detail/matrix_def.hpp"
#include "viennacl/scalar.hpp"
#include "./viennapde_enum.hpp"
// #include "viennacl/linalg/matrix_operations.hpp"
// #include "viennacl/linalg/sparse_matrix_operations.hpp"
// #include "viennacl/tools/tools.hpp"
// #include "viennacl/tools/matrix_size_deducer.hpp"
// #include "viennacl/meta/result_of.hpp"
// #include "viennacl/meta/enable_if.hpp"
// #include "viennacl/traits/handle.hpp"
// #include "viennacl/traits/row_major.hpp"






// SECTION 01a Predeclare the variable mesh class
namespace viennapde
{
template <typename NumericT>
class varmesh;

} //namespace viennapde








// SECTION 02 COPY interface with other classes
// The interface depends on the varmesh, so we need a predefinition of varmesh class declared above.
namespace viennacl
{
// SECTION 02_001 COPY interface from <- other classes
template <typename NumericT>
/** @brief Conversion: std::vector -> varmesh
 * 
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} i_std_varmesh : 
 * @param  {viennapde::varmesh<NumericT>} o_varmesh             : 
 */
void copy(  const std::vector<std::vector<std::vector<NumericT>>> & i_std_varmesh, 
            viennapde::varmesh<NumericT> *o_varmesh)
{
    size_t l_var_index_num = i_std_varmesh.size();
    for (size_t var_index = 0; var_index < l_var_index_num; var_index++)
        viennacl::copy(i_std_varmesh[var_index], o_varmesh->data_[var_index]);
}

// SECTION 02_002 COPY interface to -> other classes
template <typename NumericT>
/** @brief Conversion: varmesh -> std::vector
 * 
 * @param  {viennapde::varmesh<NumericT>} i_varmesh              : 
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} i_std_varmesh : 
 */
void copy(  const viennapde::varmesh<NumericT> & i_varmesh, 
            std::vector<std::vector<std::vector<NumericT>>> *o_std_varmesh)
{
    for (size_t var_index = 0; var_index < i_varmesh.get_var_num(); var_index++)
        viennacl::copy(i_varmesh.data_[var_index], o_std_varmesh->at(var_index));
}
} // namespace viennacl




// SECTION 01 Define the variable mesh (varmesh) class
namespace viennapde
{

template <typename NumericT>
class varmesh
{
public:
    std::vector<viennacl::matrix<NumericT>> data_;    /** @brief The data_ stored in a ViennaCL matrix organized by a STL vector */
public:
    inline size_t get_var_num()  const { return data_.size();};
    inline size_t get_row_num()    const { return data_[0].size1();};
    inline size_t get_column_num() const { return data_[0].size2();};

public:
    // SECTION 01_001 Constructor & Destructor
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} l_var_num   : 
     * @param  {ssize_t} l_row_num            : 
     * @param  {ssize_t} l_column_num : 
     */
    explicit varmesh(size_t l_var_num, ssize_t l_row_num, ssize_t l_column_num);
    /** @brief varmesh constructor 
     * @param  {std::vector<std::vector<std::vector<NumericT>>>} i_std_varmesh : 
     */
    explicit varmesh(const std::vector<std::vector<std::vector<NumericT>>> & i_std_varmesh);
    explicit varmesh(const varmesh<NumericT> & i_varmesh);
    // TODO: Move assignment not yet implemented. 
    // varmesh(varmesh<NumericT> && i_varmesh);
    /** @brief ~varmesh Destructor*/
    // ~varmesh();

    // SECTION 02_001 COPY interface from other classes
    friend void viennacl::copy<NumericT>(const std::vector<std::vector<std::vector<NumericT>>> & i_std_varmesh, viennapde::varmesh<NumericT> *o_varmesh);
    friend void viennacl::copy<NumericT>(  const viennapde::varmesh<NumericT> & i_varmesh, 
            std::vector<std::vector<std::vector<NumericT>>> *o_std_varmesh);
    
    
};

// SECTION 01_001 Constructor & Destructor
template <typename NumericT>
varmesh<NumericT>::varmesh(size_t l_var_num, ssize_t l_row_num, ssize_t l_column_num)
{
    this->data_.resize(l_var_num);
    for (size_t var_index = 0; var_index < l_var_num; var_index++)
        this->data_[var_index].resize(l_row_num, l_column_num);
}

template <typename NumericT>
varmesh<NumericT>::varmesh(const std::vector<std::vector<std::vector<NumericT>>> & i_std_varmesh)
{
    this->data_.resize(i_std_varmesh.size());
    for (size_t var_index = 0; var_index < i_std_varmesh.size(); var_index++)
    {
        this->data_[var_index].resize(i_std_varmesh[0].size(), i_std_varmesh[0][0].size());
        viennacl::copy(i_std_varmesh[var_index], this->data_[var_index]);
    }
}

template <typename NumericT>
varmesh<NumericT>::varmesh(const varmesh<NumericT> & i_varmesh)
{
    this->data_.resize(i_varmesh.get_var_num());
    for (size_t var_index = 0; var_index < i_varmesh.get_var_num(); var_index++)
    {
        this->data_[var_index].resize(i_varmesh.get_row_num(), i_varmesh.get_column_num());
        this->data_[var_index] = i_varmesh.data_[var_index];
    }
}


// template <typename NumericT>
// varmesh<NumericT>::~varmesh()
// {   
//     // delete this->data_;
// }


// SECTION 03_002a Vermesh Convolution, utilizing 01_002a
/** @brief Convolve the varmesh data by the 2D matrix kernel, which would be the base of varmesh shift and filter
 * 
 * @param  {viennacl::matrix<NumericT>} i_kernel    : 
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
            viennapde::ConvolutionType ConvolType = EQUIV,
            bool KerElementIdentity = false, 
            OptimizeLevel optimize_level = OptimizeLevel::First>
void convolve(
    const viennacl::matrix<NumericT> & i_matrix,
    const viennacl::matrix<NumericT> & i_kernel,
    viennacl::matrix<NumericT> & o_matrix,
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() ) // REVIEW: ROIrc_vec may be combined with KerElementIdentity to make the same element compute together
{
    // STUB 01 
    size_t l_kernel_size1 = (i_kernel.size1()-1)/2;
    size_t l_kernel_size2 = (i_kernel.size2()-1)/2;

    // NOTE Argument ROIrc_vec default empty case, all entries are filled in it here.
    if (ROIrc_vec.empty())
    {
        for (size_t i = 0; i < i_kernel.size1(); i++)
        for (size_t j = 0; j < i_kernel.size2(); j++)
            ROIrc_vec.push_back(std::make_pair<int, int>(i, j));
    }
    
    // STUB 02 Multiply the scalar and contribute to the final varmesh.

    o_matrix.clear(); // REVIEW This may be the efficiency bottleneck which is safe but slow 
    for(auto & iter: ROIrc_vec)
    {
        size_t  l_row = iter.first, 
                l_column = iter.second;
        int bias1 = l_row-l_kernel_size1; //TODO, here is int, make sure it will not cause leakage.
        int bias2 = l_column-l_kernel_size2;
        // REVIEW The following command is helpful to avoid unexpected too big kernel, however, I am not satisfactory with the if statement which may be a bottleneck in this loop.
        if ((std::abs(bias1)>=l_kernel_size1) || (std::abs(bias2)>=l_kernel_size2)) continue;
        if constexpr (ConvolType == ConvolutionType::EQUIV)
        {
            viennacl::range submat_row_range_from(std::max(bias1, 0), 
                                                    i_matrix.size1() + std::min(bias1, 0)),
                            submat_column_range_from(std::max(bias2, 0), 
                                                    i_matrix.size2() + std::min(bias2, 0)),
                            submat_row_range_to(std::max(-bias1, 0), 
                                                    i_matrix.size1() + std::min(-bias1, 0)),
                            submat_column_range_to(std::max(-bias2, 0), 
                                                    i_matrix.size2() + std::min(-bias2, 0));
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(i_matrix, submat_row_range_from, submat_column_range_from),
                submatrix_tobe_add(o_matrix, submat_row_range_to, submat_column_range_to); 

            submatrix_tobe_add += i_kernel(l_row, l_column) * submatrix_to_add;
        }
        else // TODO: Not yet implemented for other convolution type
        {
            std::cerr << "Not yet implemented for other convolution type than EQUIV." << std::endl;
            viennacl::range submat_row_range(l_row, l_row+i_matrix.size1()), 
                            submat_column_range(l_column, l_column+i_matrix.size2());
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to_add(i_matrix, submat_row_range, submat_column_range); 
            o_matrix += i_kernel(l_row, l_column) * submatrix_to_add;
        }
        
    }    
} //function void viennapde::convolve

// TODO totally change the comment & plut the output initialization & temp varmesh removal
// SECTION 03_002a Varmesh Convolution, utilizing 01_002a
/** @brief Convolve the varmesh data by the 2D matrix kernel, which would be the base of varmesh shift and filter
 * 
 * @param  {viennacl::matrix<NumericT>} i_kernel    : 
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
            viennapde::ConvolutionType ConvolType = EQUIV,
            bool KerElementIdentity = false, 
            OptimizeLevel optimize_level = OptimizeLevel::First>
void convolve(
    const viennapde::varmesh<NumericT> & i_varmesh,
    const viennacl::matrix<NumericT> & i_kernel,
    viennapde::varmesh<NumericT> & o_varmesh,
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() )
{
    o_varmesh.data_.resize(i_varmesh.get_var_num()); // REVIEW This may be the efficiency bottleneck which is safe but slow 
    // STUB 02 Multiply the scalar and contribute to the final varmesh.

    // NOTE Argument ROIrc_vec default empty case, all entries are filled in it here.
    if (ROIrc_vec.empty())
    {
        for (size_t i = 0; i < i_kernel.size1(); i++)
        for (size_t j = 0; j < i_kernel.size2(); j++)
            ROIrc_vec.push_back(std::make_pair<int, int>(i, j));
    }
    for (size_t var_index=0; var_index< i_varmesh.get_var_num(); var_index++) 
        viennapde::convolve<NumericT, ConvolType, KerElementIdentity, optimize_level> 
            (i_varmesh.data_[var_index], i_kernel, o_varmesh.data_[var_index], ROIrc_vec);
} //function void viennapde::convolve



// SECTION 01_002b Mesh Convolution
/** @brief the same as the convolve funciton 01_002a, just the default ROI becomes the whole matrix.
 * @param  {viennacl::matrix<NumericT>} i_kernel : 
 * 
 * @example
 * vcl_MatrixT kernel(5, 5);
 * kernel.clear();
 * kernel(2, 1) = -1.0/2; kernel(2, 2) = 0; kernel(2, 3) = 1.0/2;
 * vcl_varmesh.convolve(kernel);
 */
template <  typename NumericT, 
            viennapde::ConvolutionType ConvolType = EQUIV,
            bool KerElementIdentity = false, 
            OptimizeLevel optimize_level = OptimizeLevel::First>
void convolve(
    viennapde::varmesh<NumericT> & i_varmesh,
    const viennacl::matrix<NumericT> & i_kernel,
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() )
{
    // FIXME A strange bug here that if you use make_pair<size_t, size_t> to specify the range of interest, the compiler fails.
    viennapde::varmesh<NumericT> t_varmesh(i_varmesh);
    viennapde::convolve<NumericT, ConvolType, KerElementIdentity, optimize_level>(t_varmesh, i_kernel, i_varmesh, ROIrc_vec);
} //function varmesh::convolve




} //namespace viennapde


