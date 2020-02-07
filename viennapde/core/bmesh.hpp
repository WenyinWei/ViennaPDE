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

/** @file viennapde/core/bmesh.hpp
    @brief Varmesh with margin.
*/

#include "./mesh.hpp"




// SECTION 01 Define the variable mesh (Varmesh) class
namespace viennapde
{

template <typename NumericT>
class BVarmesh : public Varmesh<NumericT>
{
private: 
    cord3<GridIntT> margin_;
public:
    // NOTE Extending boundary will only be used while initializing, so the performance is not quite important.
    void extend_boundary(GridIntT iX, GridIntT iY, GridIntT iZ) 
    {
        // Only when the original mesh does not boundary, can the mesh extend its boundary. Otherwise we don't have corresponding functions to extend the boundary (and we don't even have plan to prepare for that). 
        assert(margin_.x==0 && margin_.y==0 && margin_.z==0);
        assert(iX >= 0 && iY >= 0 && iZ >= 0);
        margin_.x=iX; margin_.y=iY; margin_.z=iZ;
        const GridIntT  Nx = this->get_row_num(),
                        Ny = this->get_column_num(),
                        Nz = this->get_layer_num();
        const size_t tMatrixNum = 4;
        viennacl::matrix<NumericT> tMatrix[tMatrixNum]; // We prepare the temp matrixs to help facilitate the boundary extension. The number is set casually and need to be polished manually.
        viennacl::range submat_row_range_to(iX, Nx-iX), 
                        submat_column_range_to(iY, Ny-iY);        
        
        for (size_t i = 0; i < tMatrixNum; i++)  tMatrix[i].resize(Nx+2*iX, Ny+2*iY);

        for (GridIntT i = 0; i < Nz; i++)
        {
            viennacl::matrix_range<viennacl::matrix<NumericT>>  
                submatrix_to(tMatrix[i % tMatrixNum], submat_row_range_to, submat_column_range_to); 
            submatrix_to = *(this->at(i));            
            *(this->at(i)) = tMatrix[i % tMatrixNum];
        }
        
        for (GridIntT i = 0; i < iZ; i++)
        {
            this->push_front(std::make_shared<viennacl::matrix<NumericT>>((size_t)(Nx+2*iX), (size_t)(Ny+2*iY) ));
            this->push_back(std::make_shared<viennacl::matrix<NumericT>>((size_t)(Nx+2*iX), (size_t)(Ny+2*iY) ));
        }
    };
    cord3<GridIntT> get_size_boundary() const { return margin_; };
    GridIntT get_row_boundary() const { return margin_.x; };
    GridIntT get_column_boundary() const { return margin_.y; };
    GridIntT get_layer_boundary() const { return margin_.z; };

    // CTOR
    explicit BVarmesh(const DequeMat<NumericT> & iDequeMat): Varmesh<NumericT>{iDequeMat}, margin_{0,0,0} {};
    virtual ~BVarmesh() {};
    // NOTE BC functions (boundary condition) will be used frequently, so their frequency are pretty important.
    void BCPeriodic() 
    {
        const GridIntT  bx = this->get_size_boundary().x,
                        by = this->get_size_boundary().y,
                        bz = this->get_size_boundary().z;
        const GridIntT  Nx = this->get_row_num(),
                        Ny = this->get_column_num(),
                        Nz = this->get_layer_num();

        for (GridIntT i = bz; i < this->get_layer_num() - bz; i++)
        {
            viennacl::matrix<NumericT> & matThisLayer = *(this->at(i));
            {// Row Range Override
                viennacl::range submat_row_range_from(Nx-2*bx, Nx-bx),
                                submat_column_range_from(by, Ny-by),
                                submat_row_range_to(0, bx), 
                                submat_column_range_to(by, Ny-by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(matThisLayer, submat_row_range_from, submat_column_range_from),
                    submatrix_to(matThisLayer, submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
            {
                viennacl::range submat_row_range_from(bx, 2*bx),
                                submat_column_range_from(by, Ny-by),
                                submat_row_range_to(Nx-bx, Nx), 
                                submat_column_range_to(by, Ny-by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(matThisLayer, submat_row_range_from, submat_column_range_from),
                    submatrix_to(matThisLayer, submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
            {// Column Range Override
                viennacl::range submat_row_range_from(0, Nx),
                                submat_column_range_from(by, 2*by),
                                submat_row_range_to(0, Nx), 
                                submat_column_range_to(Ny-by, Ny);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(matThisLayer, submat_row_range_from, submat_column_range_from),
                    submatrix_to(matThisLayer, submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }            
            {
                viennacl::range submat_row_range_from(0, Nx),
                                submat_column_range_from(Ny-2*by, Ny-by),
                                submat_row_range_to(0, Nx), 
                                submat_column_range_to(0, by);
                viennacl::matrix_range<viennacl::matrix<NumericT>>  
                    submatrix_from(matThisLayer, submat_row_range_from, submat_column_range_from),
                    submatrix_to(matThisLayer, submat_row_range_to, submat_column_range_to); 
                submatrix_to = submatrix_from;
            }
        }

        for (GridIntT i = 0; i < bz; i++)
            *(this->at(i)) = *(this->at(i+Nz-2*bz));
        for (GridIntT i = Nz-bz; i < Nz; i++)
            *(this->at(i)) = *(this->at(i-(Nz-2*bz)));
    };

    void BCDirichlet(); //TODO Not yet implemented
    void BCNeumann(); //TODO Not yet implemented

};







} //namespace viennapde


