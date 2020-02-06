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
    using GridIntT = ssize_t;
    enum ClrOut : bool { YES=true, NO=false }; 
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
    oVarmesh->resize_ptr(iVarmeshSTL.size());
    for (size_t layer_i = 0; layer_i < iVarmeshSTL.size(); layer_i++)
    {   oVarmesh->at(layer_i).resize(iVarmeshSTL[0].size(), iVarmeshSTL[0][0].size());
        viennacl::copy(*(iVarmeshSTL[layer_i]), *(oVarmesh[layer_i]));
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
    oVarmeshSTL.resize(iVarmesh.size());
    for (size_t layer_i = 0; layer_i < iVarmesh.get_layer_num(); layer_i++)
    {
        oVarmeshSTL[layer_i].resize(iVarmesh.get_row_num());
        for (size_t row_i = 0; row_i < iVarmesh.get_row_num(); row_i++)
            oVarmeshSTL[layer_i].at(row_i).resize(iVarmesh.get_column_num());
        viennacl::copy(*(iVarmesh[layer_i]), oVarmeshSTL[layer_i]);
    }
}
} // namespace viennacl




// SECTION 01 Define the variable mesh (Varmesh) class
namespace viennapde
{
template<typename NumericT> 
using DequeMat = std::deque<std::shared_ptr<viennacl::matrix<NumericT>>>;

template <typename NumericT>
class Varmesh : public DequeMat<NumericT>
{
public:
    cord3<GridIntT> get_size_num() const { return cord3( (*this)[0]->size1(), (*this)[0]->size2(), this->size() ); };
    GridIntT get_row_num()    const { return this->at(0)->size1();};
    GridIntT get_column_num() const { return this->at(0)->size2();};
    GridIntT get_layer_num()  const { return this->size();};
public:
    /*===== SECTION Constructor & Destructor ==================================================== */
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} layer_num   : 
     * @param  {size_t} row_num     : 
     * @param  {size_t} column_num : 
     */
    explicit Varmesh(GridIntT row_num, GridIntT column_num, GridIntT layer_num): 
        DequeMat<NumericT> {(size_t)layer_num}
        {
            for (GridIntT layer_i = 0; layer_i < layer_num; layer_i++)
                this->at(layer_i) = std::make_shared<viennacl::matrix<NumericT>>(row_num, column_num);
            std::cout << "Matrix size: "<< row_num<<", "<<column_num<<"\n";
        };
    /** @brief Constructor <- VarmeshSTL
     * @param  {std::vector<std::vector<std::vector<NumericT>>>} iVarmeshSTL : 
     */
    explicit Varmesh(const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL): 
        Varmesh(iVarmeshSTL[0].size(), iVarmeshSTL[0][0].size(), iVarmeshSTL.size())
        {   
            for (GridIntT layer_i = 0; layer_i < iVarmeshSTL.size(); layer_i++)
                viennacl::copy(iVarmeshSTL[layer_i], *(this->at(layer_i)));
        };
    /** @brief Copy Constructor 
     * @param  {Varmesh<NumericT>} iDequeMat : 
     */
    explicit Varmesh(const DequeMat<NumericT> & iDequeMat, bool ref=false): 
        Varmesh(iDequeMat[0]->size1(), iDequeMat[0]->size2(), iDequeMat.size())
        {   
            for (GridIntT layer_i = 0; layer_i < iDequeMat.size(); layer_i++)
            {
                if (!ref) {
                    *(this->at(layer_i)) = *(iDequeMat[layer_i]);
                } else {
                    this->at(layer_i) = iDequeMat[layer_i];
                }
            }
        };
    explicit Varmesh(const Varmesh<NumericT> & iVarmesh, bool ref=false):Varmesh((DequeMat<NumericT>)iVarmesh) {};

    virtual ~Varmesh() {}
    
    /*===== SECTION Assignment Overriding =========================================================== */

    /** @brief Copy Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> & iVarmesh)
    {
        std::cerr << "You are using Varmesh's copy assignment which is forbidden by its unique ptr.\n" 
        << "Please check your code and debug";
        return *this;
    };
    /** @brief Move Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> && iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
        {
            this->at(i) = iVarmesh[i];
            iVarmesh[i].reset();
        }
        return *this; 
    };

    /*===== SECTION Operator Overloading=========================================================== */
    Varmesh<NumericT>& operator+ (const Varmesh<NumericT> & iVarmesh) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) += *(iVarmesh[i]);
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator- (const Varmesh<NumericT> & iVarmesh) const 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) -= *(iVarmesh[i]);
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator* (const Varmesh<NumericT> & iVarmesh) const 
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT>{this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) = viennacl::linalg::element_prod(*(this->at(i)), *(iVarmesh[i]));
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator/ (const Varmesh<NumericT> & iVarmesh) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT>{this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) = viennacl::linalg::element_div(*(this->at(i)), *(iVarmesh[i]));
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator+= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += *(iVarmesh[i]);
        return *this;
    }
    Varmesh<NumericT>& operator-= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i))  -= *(iVarmesh[i]);
        return *this;
    }
    Varmesh<NumericT>& operator*= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) = viennacl::linalg::element_prod( *(this->at(i)), *(iVarmesh[i]));
        return *this;
    }
    Varmesh<NumericT>& operator/= (const Varmesh<NumericT> & iVarmesh) 
    {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) = viennacl::linalg::element_div( *(this->at(i)), *(iVarmesh[i]));
        return *this;
    }
    Varmesh<NumericT>& operator+ (NumericT iNum) const
    {
        Varmesh<NumericT> *tVarmesh{*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) += iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator- (NumericT iNum) const
    {
        Varmesh<NumericT> *tVarmesh{*this};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) -= iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator* (NumericT iNum) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) *= iNum;
        return *tVarmesh;
    }
    Varmesh<NumericT>& operator/ (NumericT iNum) const
    {
        Varmesh<NumericT> *tVarmesh = new Varmesh<NumericT> {this->get_row_num(),this->get_column_num(),this->get_layer_num()};
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(tVarmesh->at(i)) /= iNum;
        return *tVarmesh;
    }

    /*===== SECTION COPY interface from other classes ======================================== */
    friend void viennacl::copy<NumericT>(
        const std::vector<std::vector<std::vector<NumericT>>> & iVarmeshSTL, 
        viennapde::Varmesh<NumericT> & oVarmesh);
    friend void viennacl::copy<NumericT>(
        const viennapde::Varmesh<NumericT> & iVarmesh, 
        std::vector<std::vector<std::vector<NumericT>>> & oVarmeshSTL);
    
    
    /*===== SECTION Novel Methods =========================================================== */

    /** @brief Resize the container class mesh.
     * @param  {GridIntT} Nx : 
     * @param  {GridIntT} Ny : 
     * @param  {GridIntT} Nz : 
     */
    void resize_ptr(GridIntT Nz)
    {
        GridIntT old_Nz = this->get_layer_num();
        if (Nz > old_Nz)
        {
            this->resize(Nz);
            for (GridIntT i = old_Nz; i < Nz; i++)
                this->at(i) = std::make_shared<viennacl::matrix<NumericT>>(this->get_row_num(), this->get_column_num());
        } else if (Nz = old_Nz) {
        } else if (Nz < old_Nz) {
            this->resize(Nz);
        }
    }
    
};





} //namespace viennapde


