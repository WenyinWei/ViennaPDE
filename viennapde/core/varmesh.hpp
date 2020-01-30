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

/** @file viennapde/core/Varmesh.hpp
    @brief Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
*/

// #include "viennacl/forwards.h"
// #include "viennacl/detail/matrix_def.hpp"
// #include "viennacl/scalar.hpp"
#include "./viennapde_enum.hpp"
// #include "viennacl/linalg/matrix_operations.hpp"
// #include "viennacl/linalg/sparse_matrix_operations.hpp"
// #include "viennacl/tools/tools.hpp"
// #include "viennacl/tools/matrix_size_deducer.hpp"
// #include "viennacl/meta/result_of.hpp"
// #include "viennacl/meta/enable_if.hpp"
// #include "viennacl/traits/handle.hpp"
// #include "viennacl/traits/row_major.hpp"

template <typename IntegerT>
struct cord2
{
    IntegerT x;
    IntegerT y;
    explicit cord2(IntegerT iX, IntegerT iY): x{iX}, y{iY} { };
};

template <typename IntegerT>
struct cord3
{
    IntegerT x;
    IntegerT y;
    IntegerT z = 0;
    explicit cord3(IntegerT iX, IntegerT iY, IntegerT iZ): x{iX}, y{iY}, z{iZ} { };
    explicit cord3(IntegerT iX, IntegerT iY): cord3{iX, iY, 0} { };
};



// SECTION 01a Predeclare the Variable Mesh (Varmesh) class
namespace viennapde
{
template <typename NumericT> class Varmesh;
} //namespace viennapde



// SECTION 02 COPY interface with other classes
// The interface depends on the Varmesh, so we need a predefinition of Varmesh class declared above.
namespace viennacl
{
// SECTION 02_001 COPY interface from <- other classes

/** @brief Conversion: std::vector STL -> Varmesh
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} iVarmeshSTL : 
 * @param  {viennapde::Varmesh<NumericT>} oVarmesh             : 
 */
template <typename NumericT>
void copy(  const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL, 
            viennapde::Varmesh<NumericT> & oVarmesh)
{   
    oVarmesh.data_->resize(iVarmeshSTL.size());
    for (size_t layer_i = 0; layer_i < iVarmeshSTL.size(); layer_i++)
    {   oVarmesh.data_->at(layer_i).resize(iVarmeshSTL[0].size(), iVarmeshSTL[0][0].size());
        viennacl::copy(iVarmeshSTL[layer_i], oVarmesh.data_->at(layer_i));
    }
}

// SECTION 02_002 COPY interface to -> other classes

/** @brief Conversion: Varmesh -> std::vector STL
 * @param  {viennapde::Varmesh<NumericT>} iVarmesh              : 
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} oVarmeshSTL : 
 */
template <typename NumericT>
void copy(  const viennapde::Varmesh<NumericT> & iVarmesh, 
            std::vector<std::vector<std::vector<NumericT>>> & oVarmeshSTL)
{   
    oVarmeshSTL.resize(iVarmesh.data_->size());
    for (size_t layer_i = 0; layer_i < iVarmesh.get_layer_num(); layer_i++)
    {
        oVarmeshSTL[layer_i].resize(iVarmesh.get_row_num());
        for (size_t row_i = 0; row_i < iVarmesh.get_row_num(); row_i++)
            oVarmeshSTL[layer_i].at(row_i).resize(iVarmesh.get_column_num());
        viennacl::copy(iVarmesh.data_->at(layer_i), oVarmeshSTL[layer_i]);
    }
}
} // namespace viennacl




// SECTION 01 Define the variable mesh (Varmesh) class
namespace viennapde
{

template <typename NumericT>
class Varmesh
{
    using MarginIntT = int;
public: //TODO: Too many dependency on the public, remove them first and then make it private
    std::unique_ptr<std::vector<viennacl::matrix<NumericT>>> data_;    /** @brief The data_ stored in a ViennaCL matrix organized by a STL vector */
    cord3<MarginIntT> margin_;
public:
    size_t get_layer_num()  const { return data_->size();};
    size_t get_row_num()    const { return data_->at(0).size1();};
    size_t get_column_num() const { return data_->at(0).size2();};
    cord3<MarginIntT> get_margin() const { return margin_;};
protected:
    void set_margin(MarginIntT iX, MarginIntT iY, MarginIntT iZ) {margin_.x=iX; margin_.y=iY; margin_.z=iZ;};
public:
    // SECTION 01_001 Constructor & Destructor
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} layer_num   : 
     * @param  {size_t} row_num     : 
     * @param  {size_t} column_num : 
     */
    explicit Varmesh(size_t layer_num, size_t row_num, size_t column_num): 
        data_{new std::vector<viennacl::matrix<NumericT>> (layer_num)}, margin_{0,0,0}
        {   
            for (size_t layer_i = 0; layer_i < layer_num; layer_i++)
                this->data_->at(layer_i).resize(row_num, column_num);
        };
    /** @brief Varmesh Constructor <- VarmeshSTL
     * @param  {std::vector<std::vector<std::vector<NumericT>>>} iVarmeshSTL : 
     */
    explicit Varmesh(const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL): 
        data_{new std::vector<viennacl::matrix<NumericT>> (iVarmeshSTL.size())}, margin_{0,0,0}
        {   
            for (size_t layer_i = 0; layer_i < iVarmeshSTL.size(); layer_i++)
            {   this->data_->at(layer_i).resize(iVarmeshSTL[0].size(), iVarmeshSTL[0][0].size());
                viennacl::copy(iVarmeshSTL[layer_i], this->data_->at(layer_i));
            }
        };
    /** @brief Varmesh Copy Constructor 
     * @param  {Varmesh<NumericT>} iVarmesh : 
     */
    explicit Varmesh(const Varmesh<NumericT> & iVarmesh): 
        data_{new std::vector<viennacl::matrix<NumericT>> (iVarmesh.get_layer_num())}, margin_{0,0,0}
        {   
            for (size_t layer_i = 0; layer_i < iVarmesh.get_layer_num(); layer_i++)
            {   this->data_->at(layer_i).resize(iVarmesh.get_row_num(), iVarmesh.get_column_num());
                this->data_->at(layer_i) = iVarmesh.data_->at(layer_i);
            }
        };
    // TODO: Move assignment not yet implemented. 
    // varmesh& operator= (varmesh & iVarmesh);
    // varmesh(varmesh<NumericT> && iVarmesh);
    /** @brief ~varmesh Destructor*/
    // ~varmesh();

    // SECTION 02_001 COPY interface from other classes
    friend void viennacl::copy<NumericT>(const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL, 
                                        viennapde::Varmesh<NumericT> & oVarmesh);
    friend void viennacl::copy<NumericT>(const viennapde::Varmesh<NumericT> & iVarmesh, 
                                        std::vector<std::vector<std::vector<NumericT>>> & oVarmeshSTL);
    
    
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
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() ) // REVIEW: ROIrc_vec may be combined with KerEleId to make the same element compute together
{
    // STUB 01 
    const size_t iKernelHalf1 = (iKernel.size1()-1)/2;
    const size_t iKernelHalf2 = (iKernel.size2()-1)/2;

    // NOTE Argument ROIrc_vec default empty case, all entries are filled in it here.
    if (ROIrc_vec.empty())
    {
        for (size_t i = 0; i < iKernel.size1(); i++)
        for (size_t j = 0; j < iKernel.size2(); j++)
            ROIrc_vec.push_back(std::make_pair<int, int>(i, j));
    }
    
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.

    oMatrix.clear(); // REVIEW This may be the efficiency bottleneck which is safe but slow 
    for(auto & iter: ROIrc_vec)
    {
        size_t  l_row = iter.first, 
                l_column = iter.second;
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

// TODO totally change the comment & plut the output initialization & temp Varmesh removal
// SECTION 03_002a Varmesh Convolution, utilizing 01_002a
/** @brief Convolve the Varmesh data by the 2D matrix kernel, which would be the base of Varmesh shift and filter
 * 
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
    const viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    viennapde::Varmesh<NumericT> & oVarmesh,
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() )
{
    oVarmesh.data_.resize(iVarmesh.get_layer_num()); // REVIEW This may be the efficiency bottleneck which is safe but slow 
    // STUB 02 Multiply the scalar and contribute to the final Varmesh.

    // NOTE Argument ROIrc_vec default empty case, all entries are filled in it here.
    if (ROIrc_vec.empty())
    {
        for (size_t i = 0; i < iKernel.size1(); i++)
        for (size_t j = 0; j < iKernel.size2(); j++)
            ROIrc_vec.push_back(std::make_pair<int, int>(i, j));
    }
    for (size_t layer_i=0; layer_i< iVarmesh.get_layer_num(); layer_i++) 
        viennapde::convolve<NumericT, ConvolType> 
            (iVarmesh.data_[layer_i], iKernel, oVarmesh.data_[layer_i], ROIrc_vec);
} //function void viennapde::convolve



// SECTION 01_002b Mesh Convolution
/** @brief the same as the convolve funciton 01_002a, just the default ROI becomes the whole matrix.
 * @param  {viennacl::matrix<NumericT>} iKernel : 
 * 
 * @example
 * vcl_MatrixT kernel(5, 5);
 * kernel.clear();
 * kernel(2, 1) = -1.0/2; kernel(2, 2) = 0; kernel(2, 3) = 1.0/2;
 * vcl_varmesh.convolve(kernel);
 */
template <  typename NumericT, 
            viennapde::ConvolutionType ConvolType = EQUIV>
void convolve(
    viennapde::Varmesh<NumericT> & iVarmesh,
    const viennacl::matrix<NumericT> & iKernel,
    std::vector<std::pair<size_t, size_t>> ROIrc_vec = std::vector<std::pair<size_t, size_t>>() )
{
    // FIXME A strange bug here that if you use make_pair<size_t, size_t> to specify the range of interest, the compiler fails.
    viennapde::Varmesh<NumericT> t_varmesh(iVarmesh);
    viennapde::convolve<NumericT, ConvolType>(t_varmesh, iKernel, iVarmesh, ROIrc_vec);
} //function Varmesh::convolve




} //namespace viennapde


