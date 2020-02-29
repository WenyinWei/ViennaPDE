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

/** @file viennapde/core/mesh.hpp
    @brief Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
*/

#include <vector> // viennacl::copy Vec<Vec<Vec<>>><->DequeMat
#include <deque> // using DequeMat = std::deque<std::shared_ptr<viennacl::matrix<NumericT>>>;
#include <memory> // std::shared_ptr

// ViennaCL library
#include "viennacl/matrix.hpp"
#include "viennacl/detail/matrix_def.hpp"
// #include "viennacl/linalg/matrix_operations.hpp"
// #include "viennacl/linalg/sparse_matrix_operations.hpp"
// #include "viennacl/tools/tools.hpp"
// #include "viennacl/tools/matrix_size_deducer.hpp"
// #include "viennacl/meta/result_of.hpp"
// #include "viennacl/meta/enable_if.hpp"
// #include "viennacl/traits/handle.hpp"
// #include "viennacl/traits/row_major.hpp"
#include "./cord.hpp"



// SECTION 01a Predeclare the Variable Mesh (mesh) class
namespace viennapde
{
    template <typename NumericT> class mesh;
    template <typename NumericT> class meshb;
    using GridIntT = ssize_t;
    enum ClrOut : bool { YES=true, NO=false };
} //namespace viennapde



// SECTION 02 COPY interface with other classes
// The interface depends on the mesh, so we need a predefinition of mesh class declared above.
namespace viennacl
{
// SECTION 02_001 COPY interface from <- other classes

/** @brief Conversion: std::vector STL -> mesh
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} iMeshSTL : 
 * @param  {viennapde::mesh<NumericT>} oMesh             : 
 */
template <typename NumericT>
void copy(  const std::vector<std::vector<std::vector<NumericT>>> & iMeshSTL, 
            viennapde::mesh<NumericT> & oMesh)
{   
    oMesh.resize_ptr(iMeshSTL.size());
    for (size_t layer_i = 0; layer_i < iMeshSTL.size(); layer_i++)
    {   oMesh[layer_i]->resize(iMeshSTL[0].size(), iMeshSTL[0][0].size());
        viennacl::copy(iMeshSTL[layer_i], *(oMesh[layer_i]));
    }
}

// SECTION 02_002 COPY interface to -> other classes

/** @brief Conversion: mesh -> std::vector STL
 * @param  {viennapde::mesh<NumericT>} iMesh              : 
 * @param  {std::vector<std::vector<std::vector<NumericT>>>} oMeshSTL : 
 */
template <typename NumericT>
void copy(  const viennapde::mesh<NumericT> & iMesh, 
            std::vector<std::vector<std::vector<NumericT>>> & oMeshSTL)
{   
    oMeshSTL.resize(iMesh.size());
    for (size_t layer_i = 0; layer_i < iMesh.get_layer_num(); layer_i++)
    {
        oMeshSTL[layer_i].resize(iMesh.get_row_num());
        for (size_t row_i = 0; row_i < iMesh.get_row_num(); row_i++)
            oMeshSTL[layer_i].at(row_i).resize(iMesh.get_column_num());
        viennacl::copy(*(iMesh[layer_i]), oMeshSTL[layer_i]);
    }
}
} // namespace viennacl




// SECTION 01 Define the variable mesh (mesh) class
namespace viennapde
{
template<typename NumericT> 
using DequeMat = std::deque<std::shared_ptr<viennacl::matrix<NumericT>>>;

template <typename NumericT>
class mesh : public DequeMat<NumericT>
{
    friend class meshb<NumericT>;
public:
    cord3<size_t> get_size_num() const { return cord3(get_row_num(), get_column_num(), this->size() ); };
    size_t get_row_num()    const { return this->empty() ? 0 : this->at(0)->size1();};
    size_t get_column_num() const { return this->empty() ? 0 : this->at(0)->size2();};
    size_t get_layer_num()  const { return this->size();};
public:
    /*===== SECTION Constructor & Destructor ==================================================== */
    template <class InputIterator>
    mesh(
        InputIterator first, 
        InputIterator last) : DequeMat<NumericT>{first, last} {}; // @brief Range CTOR inherited from std::deque
    explicit mesh(cord3<size_t> size_cord3): mesh(size_cord3.x, size_cord3.y, size_cord3.z) {};
    /** @brief Constuctor for the varmesh class by an existing 3D std::vector class
     * @param  {size_t} layer_num   : 
     * @param  {size_t} row_num     : 
     * @param  {size_t} column_num : 
     */
    explicit mesh(size_t row_num, size_t column_num, size_t layer_num): 
        DequeMat<NumericT> {layer_num}
        {
            for (GridIntT layer_i = 0; layer_i < layer_num; layer_i++)
                this->at(layer_i) = std::make_shared<viennacl::matrix<NumericT>>(row_num, column_num);
        };
    /** @brief Constructor <- VarmeshSTL
     * @param  {std::vector<std::vector<std::vector<NumericT>>>} iMeshSTL : 
     */
    explicit mesh(const std::vector<std::vector<std::vector<NumericT>>> & iMeshSTL): 
        mesh(iMeshSTL[0].size(), iMeshSTL[0][0].size(), iMeshSTL.size())
        {   
            for (GridIntT layer_i = 0; layer_i < iMeshSTL.size(); layer_i++)
                viennacl::copy(iMeshSTL[layer_i], *(this->at(layer_i)));
        };
    /** @brief Copy Constructor */
    mesh(const mesh<NumericT> & iMesh):
        DequeMat<NumericT> {iMesh.size()}
        {   
            for (GridIntT layer_i = 0; layer_i < iMesh.size(); layer_i++)
            {
                this->at(layer_i) = std::make_shared<viennacl::matrix<NumericT>>();
                *(this->at(layer_i)) = *(iMesh[layer_i]);
            }
        };
    /** @brief  Move Constructor */
    mesh(mesh<NumericT> && iMesh):
        DequeMat<NumericT> {iMesh.size()}
        {   
            for (GridIntT layer_i = 0; layer_i < iMesh.size(); layer_i++)
            {
                this->at(layer_i) = iMesh[layer_i];
                iMesh[layer_i].reset();
            }
            iMesh.resize(0);
            iMesh.shrink_to_fit();
        };
    virtual ~mesh() {}
    
    /*===== SECTION Assignment Overriding =========================================================== */

    /** @brief Copy Assignment */
    mesh<NumericT>& operator= (mesh<NumericT> & iMesh) = delete; // Only allow copy constructor to work for copying
    /** @brief Move Assignment */
    mesh<NumericT>& operator= (mesh<NumericT> && iMesh) {
        this->resize(this->get_layer_num());
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            this->at(i) = iMesh[i];
        iMesh.resize(0);
        iMesh.shrink_to_fit();
        return *this; 
    };
    mesh<NumericT>& operator+= (const mesh<NumericT> & iMesh) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += *(iMesh[i]);
        return *this;
    }
    mesh<NumericT>& operator-= (const mesh<NumericT> & iMesh) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i))  -= *(iMesh[i]);
        return *this;
    }
    mesh<NumericT>& operator*= (const mesh<NumericT> & iMesh) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) = viennacl::linalg::element_prod( *(this->at(i)), *(iMesh[i]));
        return *this;
    }
    mesh<NumericT>& operator/= (const mesh<NumericT> & iMesh) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) = viennacl::linalg::element_div( *(this->at(i)), *(iMesh[i]));
        return *this;
    }
    mesh<NumericT>& operator+= (const char & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator-= (const char & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) -= scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator*= (const char & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) *= iNum;
        return *this;
    }
    mesh<NumericT>& operator/= (const char & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) /= iNum;
        return *this;
    }
    mesh<NumericT>& operator+= (const short & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator-= (const short & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) -= scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator*= (const short & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) *= iNum;
        return *this;
    }
    mesh<NumericT>& operator/= (const short & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) /= iNum;
        return *this;
    }
    mesh<NumericT>& operator+= (const int & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator-= (const int & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) -= scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator*= (const int & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) *= iNum;
        return *this;
    }
    mesh<NumericT>& operator/= (const int & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) /= iNum;
        return *this;
    }
    mesh<NumericT>& operator+= (const float & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator-= (const float & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) -= scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator*= (const float & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) *= iNum;
        return *this;
    }
    mesh<NumericT>& operator/= (const float & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) /= iNum;
        return *this;
    }
    mesh<NumericT>& operator+= (const double & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) += scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator-= (const double & iNum) {
        viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(this->get_row_num(), this->get_column_num(), iNum);
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) -= scalar_mat;
        return *this;
    }
    mesh<NumericT>& operator*= (const double & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) *= iNum;
        return *this;
    }
    mesh<NumericT>& operator/= (const double & iNum) {
        for (GridIntT i = 0; i < this->get_layer_num(); i++)
            *(this->at(i)) /= iNum;
        return *this;
    }

    /*===== SECTION COPY interface from other classes ======================================== */
    friend void viennacl::copy<NumericT>(
        const std::vector<std::vector<std::vector<NumericT>>> & iMeshSTL, 
        viennapde::mesh<NumericT> & oMesh);
    friend void viennacl::copy<NumericT>(
        const viennapde::mesh<NumericT> & iMesh, 
        std::vector<std::vector<std::vector<NumericT>>> & oMeshSTL);
    
    
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
        bool is_empty = this->empty();
        if (Nz > old_Nz)
        {
            this->resize(Nz);
            if (((!is_empty) && Nx==0 && Ny==0) 
            || ((!is_empty) && Nx==this->get_row_num() && Ny==this->get_column_num()))
            {
                for (GridIntT i = old_Nz; i < Nz; i++)
                    this->at(i) = std::make_shared<viennacl::matrix<NumericT>>(this->get_row_num(), this->get_column_num());
            } else if ((!is_empty) && Nx!=0 && Ny!=0 && Nx!=this->get_row_num() && Ny!=this->get_column_num())
            {
                std::cerr << "Mesh::Resize_ptr function, original size is different from your input matrix size.\n"
                    << "Please check your mesh size or use clear method first.\n";
            } else if (is_empty && Nx!=0 && Ny!=0)
            {
                for (GridIntT i = old_Nz; i < Nz; i++)
                    this->at(i) = std::make_shared<viennacl::matrix<NumericT>>(Nx, Ny);
            } else if (is_empty && Nx==0 && Ny==0)
            {
                throw("Mesh::Resize_ptr function, please tell how big a matrix you want.");
                std::cerr << "Mesh::Resize_ptr function, please tell how big a matrix you want.\n";
            } else {std::cerr << "Mesh::Resize_ptr function, your input is not legal.\n";}
        } else  { //if (Nz <= old_Nz)
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
inline mesh<NumericT> operator+ (mesh<NumericT> lhs, const mesh<NumericT> & rhs) {
    lhs += rhs;
    return lhs;
}
template <typename NumericT>
inline mesh<NumericT> operator- (mesh<NumericT> lhs, const mesh<NumericT> & rhs) {
    lhs -= rhs;
    return lhs;
}
template <typename NumericT>
inline mesh<NumericT> operator* (mesh<NumericT> lhs, const mesh<NumericT> & rhs) {
    lhs *= rhs;
    return lhs;
}
template <typename NumericT>
inline mesh<NumericT> operator/ (mesh<NumericT> lhs, const mesh<NumericT> & rhs) {
    lhs /= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator+ (mesh<NumericT_l> lhs, NumericT_r rhs) {
    lhs += rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator- (mesh<NumericT_l> lhs, NumericT_r rhs) {
    lhs -= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator* (mesh<NumericT_l> lhs, NumericT_r rhs) {
    lhs *= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator/ (mesh<NumericT_l> lhs, NumericT_r rhs) {
    lhs /= rhs;
    return lhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator+ (NumericT_l lhs, mesh<NumericT_r> rhs) {
    rhs += lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator- (NumericT_l lhs, mesh<NumericT_r> rhs) {
    rhs -= lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator* (NumericT_l lhs, mesh<NumericT_r> rhs) {
    rhs *= lhs;
    return rhs;
}
template <typename NumericT_l, typename NumericT_r>
inline mesh<NumericT_l> operator/ (NumericT_l lhs, mesh<NumericT_r> rhs) {
    rhs /= lhs;
    return rhs;
}
} //namespace viennapde



