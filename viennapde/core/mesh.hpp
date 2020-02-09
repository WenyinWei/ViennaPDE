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

#include <vector> // viennacl::copy Vec<Vec<Vec<>>><->DequeMat
#include <deque> // using DequeMat = std::deque<std::shared_ptr<viennacl::matrix<NumericT>>>;
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
    cord3<GridIntT> get_size_num() const { return cord3( (GridIntT) this->at(0)->size1(), (GridIntT) this->at(0)->size2(), (GridIntT) this->size() ); };
    GridIntT get_row_num()    const { return this->at(0)->size1();};
    GridIntT get_column_num() const { return this->at(0)->size2();};
    GridIntT get_layer_num()  const { return this->size();};
public:
    /*===== SECTION Constructor & Destructor ==================================================== */

    explicit Varmesh(
        const typename DequeMat<NumericT>::iterator & first, 
        const typename DequeMat<NumericT>::iterator & last) : DequeMat<NumericT>{first, last} {}; // @brief Range CTOR inherited from std::deque
    explicit Varmesh(cord3<GridIntT> size_cord3): Varmesh(size_cord3.x, size_cord3.y, size_cord3.z) {};
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} layer_num   : 
     * @param  {size_t} row_num     : 
     * @param  {size_t} column_num : 
     */
    explicit Varmesh(size_t row_num, size_t column_num, size_t layer_num): 
        DequeMat<NumericT> {layer_num}
        {
            for (GridIntT layer_i = 0; layer_i < layer_num; layer_i++)
                this->at(layer_i) = std::make_shared<viennacl::matrix<NumericT>>(row_num, column_num);
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
    Varmesh(const DequeMat<NumericT> & iDequeMat):
        DequeMat<NumericT> {iDequeMat.size()}
        {   
            for (GridIntT layer_i = 0; layer_i < iDequeMat.size(); layer_i++)
            {
                this->at(layer_i) = std::make_shared<viennacl::matrix<NumericT>>();
                *(this->at(layer_i)) = *(iDequeMat[layer_i]);
            }
        };
    Varmesh(const Varmesh<NumericT> & iVarmesh): Varmesh((DequeMat<NumericT>)iVarmesh) {};
    /** @brief  Move Constructor */
    Varmesh(Varmesh<NumericT> && iVarmesh):
        DequeMat<NumericT> {iVarmesh.size()}
        {   
            for (GridIntT layer_i = 0; layer_i < iVarmesh.size(); layer_i++)
            {
                this->at(layer_i) = iVarmesh[layer_i];
                iVarmesh[layer_i].reset();
            }
            iVarmesh.resize(0);
            iVarmesh.shrink_to_fit();
        };
    virtual ~Varmesh() {}
    
    /*===== SECTION Assignment Overriding =========================================================== */

    /** @brief Copy Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> & iVarmesh) = delete; // Only allow copy constructor to work for copying
    /** @brief Move Assignment */
    Varmesh<NumericT>& operator= (Varmesh<NumericT> && iVarmesh) 
    {
        this->resize(this->get_layer_num());
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->at(i) = iVarmesh[i];
        iVarmesh.resize(0);
        iVarmesh.shrink_to_fit();
        return *this; 
    };
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
    void resize_ptr(GridIntT Nz, GridIntT Nx=0, GridIntT Ny=0)
    {
        assert(Nx>=0 && Ny>=0 && Nz>=0);
        GridIntT old_Nz = this->get_layer_num();
        if (Nz > old_Nz)
        {
            this->resize(Nz);
            if ((!(this->empty()) && Nx==0 && Ny==0) || (!(this->empty()) && Nx==this->get_row_num() && Ny==this->get_column_num()))
            {
                for (GridIntT i = old_Nz; i < Nz; i++)
                    this->at(i) = std::make_shared<viennacl::matrix<NumericT>>(this->get_row_num(), this->get_column_num());
            } else if (!(this->empty()) && Nx!=0 && Ny!=0 && Nx!=this->get_row_num() && Ny!=this->get_column_num())
            {
                std::cerr << "Mesh::Resize_ptr function, original size is different from your input matrix size.\n";
            } else if (this->empty() && Nx!=0 && Ny!=0)
            {
                for (GridIntT i = old_Nz; i < Nz; i++)
                    this->at(i) = std::make_shared<viennacl::matrix<NumericT>>(Nx, Ny);
            } else if (this->empty() && Nx==0 && Ny==0)
            {
                throw("Mesh::Resize_ptr function, please tell how big a matrix you want.");
                std::cerr << "Mesh::Resize_ptr function, please tell how big a matrix you want.\n";
            } else {std::cerr << "Mesh::Resize_ptr function, your input is not legal.\n";}
        } else if (Nz <= old_Nz) {
            this->resize(Nz); 
        }
    }
    
    void clear()
    {
        for (size_t layer_i=0; layer_i< this->get_layer_num(); layer_i++)  this->at(layer_i)->clear();
    }
};



/*===== SECTION Operator Overloading ======================================== */
template <typename NumericT>
inline Varmesh<NumericT> operator+ (Varmesh<NumericT> lhs, const Varmesh<NumericT> & rhs) 
{
    lhs += rhs;
    return lhs;
}
template <typename NumericT>
inline Varmesh<NumericT> operator- (Varmesh<NumericT> lhs, const Varmesh<NumericT> & rhs) 
{
    lhs -= rhs;
    return lhs;
}
template <typename NumericT>
inline Varmesh<NumericT> operator* (Varmesh<NumericT> lhs, const Varmesh<NumericT> & rhs) 
{
    lhs *= rhs;
    return lhs;
}
template <typename NumericT>
inline Varmesh<NumericT> operator/ (Varmesh<NumericT> lhs, const Varmesh<NumericT> & rhs) 
{
    lhs /= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator+ (Varmesh<NumericT_l> lhs, NumericT_r rhs) 
{
    for (GridIntT i = 0; i < lhs.get_layer_num(); i++)
        *(lhs[i]) += rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator- (Varmesh<NumericT_l> lhs, NumericT_r rhs) 
{
    for (GridIntT i = 0; i < lhs.get_layer_num(); i++)
        *(lhs[i]) -= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator* (Varmesh<NumericT_l> lhs, NumericT_r rhs) 
{
    for (GridIntT i = 0; i < lhs.get_layer_num(); i++)
        *(lhs[i]) *= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator/ (Varmesh<NumericT_l> lhs, NumericT_r rhs) 
{
    for (GridIntT i = 0; i < lhs.get_layer_num(); i++)
        *(lhs[i]) /= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator+ (NumericT_l lhs, Varmesh<NumericT_r> rhs) 
{
    for (GridIntT i = 0; i < rhs.get_layer_num(); i++)
        *(rhs[i]) += lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator- (NumericT_l lhs, Varmesh<NumericT_r> rhs) 
{
    for (GridIntT i = 0; i < rhs.get_layer_num(); i++)
        *(rhs[i]) -= lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator* (NumericT_l lhs, Varmesh<NumericT_r> rhs) 
{
    for (GridIntT i = 0; i < rhs.get_layer_num(); i++)
        *(rhs[i]) *= lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline Varmesh<NumericT_l> operator/ (NumericT_l lhs, Varmesh<NumericT_r> rhs) 
{
    for (GridIntT i = 0; i < rhs.get_layer_num(); i++)
        *(rhs[i]) /= lhs;
    return rhs;
}
} //namespace viennapde


