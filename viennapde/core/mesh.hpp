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

#include <vector>
#include <deque>
#include <memory> // std::shared_ptr

// #include "viennacl/linalg/matrix_operations.hpp"
// #include "viennacl/linalg/sparse_matrix_operations.hpp"
// #include "viennacl/tools/tools.hpp"
// #include "viennacl/tools/matrix_size_deducer.hpp"
// #include "viennacl/meta/result_of.hpp"
// #include "viennacl/meta/enable_if.hpp"
// #include "viennacl/traits/handle.hpp"
// #include "viennacl/traits/row_major.hpp"
#include "./cord.hpp"



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
class Varmesh : public std::deque<std::shared_ptr<viennacl::matrix<NumericT>>> 
{
public: //TODO: Too many dependency on the public, remove them first and then make it private
    std::unique_ptr<std::deque<viennacl::matrix<NumericT>>> data_;    /** @brief The data_ stored in a ViennaCL matrix organized by a STL vector */
public:
    cord3<GridIntT> get_size_num() const { return cord3( data_->at(0).size1(), data_->at(0).size2(), data_->size() ); };
    GridIntT get_layer_num()  const { return data_->size();};
    GridIntT get_row_num()    const { return data_->at(0).size1();};
    GridIntT get_column_num() const { return data_->at(0).size2();};
public:
    // SECTION 01_001 Constructor & Destructor
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} layer_num   : 
     * @param  {size_t} row_num     : 
     * @param  {size_t} column_num : 
     */
    explicit Varmesh(GridIntT row_num, GridIntT column_num, GridIntT layer_num): 
        data_{new std::deque<viennacl::matrix<NumericT>> (layer_num)}
        {   
            for (size_t layer_i = 0; layer_i < layer_num; layer_i++)
                this->data_->at(layer_i).resize(row_num, column_num);
        };
    /** @brief Constructor <- VarmeshSTL
     * @param  {std::vector<std::vector<std::vector<NumericT>>>} iVarmeshSTL : 
     */
    explicit Varmesh(const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL): 
        Varmesh(iVarmeshSTL.size(), iVarmeshSTL[0].size(), iVarmeshSTL[0][0].size())
        {   
            for (size_t layer_i = 0; layer_i < iVarmeshSTL.size(); layer_i++)
                viennacl::copy(iVarmeshSTL[layer_i], this->data_->at(layer_i));
        };
    /** @brief Copy Constructor 
     * @param  {Varmesh<NumericT>} iVarmesh : 
     */
    explicit Varmesh(const Varmesh<NumericT> & iVarmesh, bool blankMesh=false): 
        Varmesh(iVarmesh.get_layer_num(), iVarmesh.get_row_num(), iVarmesh.get_column_num())
        {   
            if (!blankMesh)
            for (size_t layer_i = 0; layer_i < iVarmesh.get_layer_num(); layer_i++)
                this->data_->at(layer_i) = iVarmesh.data_->at(layer_i);
        };
    /** @brief Copy Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> & iVarmesh)
    {
        std::cerr << "You are using Varmesh's copy assignment which is forbidden by its unique ptr.\n" 
        << "Please check your code and debug";
        return *this;
    };
    /** Move Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> && iVarmesh) { this->data_ = std::move(iVarmesh.data_); return *this; };
    // varmesh(varmesh<NumericT> && iVarmesh);
    /** @brief ~varmesh Destructor */
    // ~varmesh();

    Varmesh<NumericT>& operator+ (const Varmesh<NumericT> & iVarmesh) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) += iVarmesh.data_->at(i);
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator- (const Varmesh<NumericT> & iVarmesh) const 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) -= iVarmesh.data_->at(i);
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator* (const Varmesh<NumericT> & iVarmesh) const 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) = viennacl::linalg::element_prod(this->data_->at(i), iVarmesh.data_->at(i));
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator/ (const Varmesh<NumericT> & iVarmesh) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) = viennacl::linalg::element_div(this->data_->at(i), iVarmesh.data_->at(i));
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator+= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->data_->at(i) += iVarmesh.data_->at(i);
        return *this;
    }
    Varmesh<NumericT>& operator-= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->data_->at(i) -= iVarmesh.data_->at(i);
        return *this;
    }
    Varmesh<NumericT>& operator*= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->data_->at(i) = viennacl::linalg::element_prod( this->data_->at(i), iVarmesh.data_->at(i));
        return *this;
    }
    Varmesh<NumericT>& operator/= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->data_->at(i) = viennacl::linalg::element_div( this->data_->at(i), iVarmesh.data_->at(i));
        return *this;
    }
    Varmesh<NumericT>& operator+ (NumericT iNum) 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) += iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator- (NumericT iNum) 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) -= iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator* (NumericT iNum) 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) *= iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator/ (NumericT iNum) 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            tVarmesh->data_->at(i) /= iNum;
        return *tVarmesh;
    }

    // SECTION 02_001 COPY interface from other classes
    friend void viennacl::copy<NumericT>(
        const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL, 
        viennapde::Varmesh<NumericT> & oVarmesh);
    friend void viennacl::copy<NumericT>(
        const viennapde::Varmesh<NumericT> & iVarmesh, 
        std::vector<std::vector<std::vector<NumericT>>> & oVarmeshSTL);
    
    
};





} //namespace viennapde


