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

/** @file viennapde/core/meshb.hpp
    @brief mesh with margin.
*/

#include "viennapde/core/mesh.hpp"




// SECTION 01 Define the variable mesh (mesh) class
namespace viennapde
{

enum BoundaryCondition{
    Unknown,
    Periodic,
    Neumann,
    Dirichlet,
};

template <typename NumericT>
class meshb 
{
private: 
    viennapde::mesh<NumericT> & mesh_ref_; 
    cord3<size_t> margin_;
    // cord3<BoundaryCondition> BC_;
    // std::weak_ptr<meshb> x_neg, x_pos, y_neg, y_pos, z_neg, z_pos;
public:
    meshb(mesh<NumericT> & iMesh): mesh_ref_{iMesh}, margin_{0,0,0} {iMesh.bdry_ptr_ = this;}; //@brief CTOR
    ~meshb() {mesh_ref_.bdry_ptr_ = nullptr;}; // @ DTOR

    cord3<size_t> get_margin() {return cord3(get_marginx(), get_marginy(), get_marginz());}
    size_t get_marginx() {return margin_.x;}
    size_t get_marginy() {return margin_.y;}
    size_t get_marginz() {return margin_.z;}

    // NOTE Extending boundary will only be used while initializing, so the performance is not quite important.
    void extend(size_t iX, size_t iY, size_t iZ) 
    {
        // Only when the original mesh does not have boundary, can the mesh extend its boundary. Otherwise we don't have corresponding functions to extend the boundary (and that's not in future plan). 
        assert(margin_.x==0 && margin_.y==0 && margin_.z==0);
        assert(iX >= 0 && iY >= 0 && iZ >= 0);
        margin_.x=iX; margin_.y=iY; margin_.z=iZ;
        const size_t Nx = mesh_ref_.get_row_num(),
                     Ny = mesh_ref_.get_column_num(),
                     Nz = mesh_ref_.get_layer_num();
        const size_t tMatrixNum = 4;
        viennacl::matrix<NumericT> tMatrix[tMatrixNum]; // We prepare the temp matrixs to help facilitate the boundary extension. The number is set casually and need to be polished manually.
        viennacl::range submat_row_range_to(iX, Nx+iX), 
                        submat_column_range_to(iY, Ny+iY);        
        
        for (size_t i = 0; i < tMatrixNum; i++)  tMatrix[i].resize(Nx+2*iX, Ny+2*iY);

        for (size_t i = 0; i < Nz; i++)
        {
            mesh_ref_[i]->resize(Nx+2*iX, Ny+2*iY);
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to(tMatrix[i % tMatrixNum], submat_row_range_to, submat_column_range_to); 
            submatrix_to = *(mesh_ref_[i]);            
            *(mesh_ref_[i]) = tMatrix[i % tMatrixNum];
        }
        
        for (size_t i = 0; i < iZ; i++)
        {
            // TODO: Check the following command is faster.
            // mesh_ref_.emplace_front((size_t)(Nx+2*iX), (size_t)(Ny+2*iY));
            // mesh_ref_.emplace_back((size_t)(Nx+2*iX), (size_t)(Ny+2*iY));
            mesh_ref_.push_front(std::make_shared<viennacl::matrix<NumericT>>((size_t)(Nx+2*iX), (size_t)(Ny+2*iY) ));
            mesh_ref_.push_back(std::make_shared<viennacl::matrix<NumericT>>((size_t)(Nx+2*iX), (size_t)(Ny+2*iY) ));
        }
    };

    // NOTE BC functions (boundary condition) will be used frequently, so their performance are pretty important.
    void refresh() 
    {
        const size_t bx = margin_.x,
                     by = margin_.y,
                     bz = margin_.z;
        const size_t Nx = mesh_ref_.get_row_num(),
                     Ny = mesh_ref_.get_column_num(),
                     Nz = mesh_ref_.get_layer_num();

        for (size_t i = bz; i < Nz - bz; i++)
        {
            {// Row Boundary Data Transfer
                viennacl::range submat_row_range_from(Nx-2*bx, Nx-bx),
                                submat_column_range_from(by, Ny-by),
                                submat_row_range_to(0, bx), 
                                submat_column_range_to(by, Ny-by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(*(mesh_ref_[i]), submat_row_range_from, submat_column_range_from),
                    submatrix_to(*(mesh_ref_[i]), submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
            {
                viennacl::range submat_row_range_from(bx, 2*bx),
                                submat_column_range_from(by, Ny-by),
                                submat_row_range_to(Nx-bx, Nx), 
                                submat_column_range_to(by, Ny-by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(*(mesh_ref_[i]), submat_row_range_from, submat_column_range_from),
                    submatrix_to(*(mesh_ref_[i]), submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
            {// Column Boundary Data Transfer
                viennacl::range submat_row_range_from(0, Nx),
                                submat_column_range_from(by, 2*by),
                                submat_row_range_to(0, Nx), 
                                submat_column_range_to(Ny-by, Ny);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(*(mesh_ref_[i]), submat_row_range_from, submat_column_range_from),
                    submatrix_to(*(mesh_ref_[i]), submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }            
            {
                viennacl::range submat_row_range_from(0, Nx),
                                submat_column_range_from(Ny-2*by, Ny-by),
                                submat_row_range_to(0, Nx), 
                                submat_column_range_to(0, by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(*(mesh_ref_[i]), submat_row_range_from, submat_column_range_from),
                    submatrix_to(*(mesh_ref_[i]), submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
        }
        // Layer Boundary Data Transfer
        for (size_t i = 0; i < bz; i++)
            *(mesh_ref_[i]) = *(mesh_ref_[i+Nz-2*bz]);
        for (size_t i = Nz-bz; i < Nz; i++)
            *(mesh_ref_[i]) = *(mesh_ref_[i-(Nz-2*bz)]);
    };



};


} //namespace viennapde


