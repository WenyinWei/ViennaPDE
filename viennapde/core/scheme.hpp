#pragma once
/* =========================================================================
   Copyright (c) 2016-2020, Department of Engineering Physics,
                            Tsinghua University, Beijing, China.

   Portions of this software are copyright by UChicago Argonne, LLC and ViennaCL team.

                            -----------------
                  ViennaPDE - The Vienna PDE Library
                            -----------------

   Project Head:    Wenyin Wei                   weiwy16@mails.tsinghua.edu.cn

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennapde/core/scheme.hpp
    @brief Implementation of concrete scheme for PDE.
*/
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/bmesh.hpp"
#include "./convol_mesh.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

template <typename NumericT>
void Godunov(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    viennapde::Varmesh<NumericT> & oVarmesh, NumericT dt, NumericT dx, unsigned int shift=1)
{
    std::cout << "I am fine 1\n";
    viennapde::BVarmesh<NumericT> tVarmesh{iVarmesh};
    std::cout << "I am fine 2\n";
    tVarmesh.extend_boundary(1,0,0);
    std::cout << "I am fine 3\n";
    tVarmesh.BCPeriodic();
    std::cout << "I am fine 4\n";
    // ssize_t shift = 1;
    // viennapde::Varmesh<NumericT> tVarmesh(1, iVarmesh.get_row_num(), iVarmesh.get_column_num()+2*shift);
    // viennacl::range submat_row_range_from(0, 1),
    //                 submat_column_range_from(0, iVarmesh.get_column_num()),
    //                 submat_row_range_to(0, 1),
    //                 submat_column_range_to(shift, iVarmesh.get_column_num()+shift);
    // viennacl::matrix_range<viennacl::matrix<NumericT>>  
    //     submatrix_to_replace(iVarmesh.data_->at(0), submat_row_range_from, submat_column_range_from),
    //     submatrix_tobe_replace(tVarmesh.data_->at(0), submat_row_range_to, submat_column_range_to); 
    // submatrix_tobe_replace = submatrix_to_replace;
    // for (ssize_t i = 0; i < shift; i++)
    // {
    //     tVarmesh.data_->at(0)(0, i) = iVarmesh.data_->at(0)(0, iVarmesh.get_column_num()+i);
    //     tVarmesh.data_->at(0)(0, iVarmesh.get_column_num()+2*shift-1-i) = iVarmesh.data_->at(0)(0, 2*shift-1-i);
    // }
    


    // Calculate the u_{j+1/2}^-, u_{j+1/2}^+
    std::vector<std::vector<std::vector<NumericT>>> t_std_varmesh_l, t_std_varmesh_r,
                                                    t_std_ustar_varmesh;
    std::cout << "I am fine 5\n";
    
    viennacl::copy(tVarmesh, t_std_varmesh_l);
    viennacl::copy(tVarmesh, t_std_varmesh_r);
    viennacl::copy(tVarmesh, t_std_ustar_varmesh);

    std::cout << "I am fine 6\n";

    // Calculate the u^*
    for (size_t i = 0; i < tVarmesh.get_layer_num(); i++)
    for (size_t j = 0; j < tVarmesh.get_row_num(); j++)
    {
            for (ssize_t k = 0; k < tVarmesh.get_column_num()-1; k++)
    {
        if ((t_std_varmesh_l[i][j][k]>=0) && (t_std_varmesh_r[i][j][k+1]>=0)) 
        {
            t_std_ustar_varmesh[i][j][k]=t_std_varmesh_l[i][j][k];
        }
        else if ((t_std_varmesh_l[i][j][k]<=0) && (t_std_varmesh_r[i][j][k+1]<=0))
        {
            t_std_ustar_varmesh[i][j][k]=t_std_varmesh_r[i][j][k+1];
        }
        else if ((t_std_varmesh_l[i][j][k]>=0) && (t_std_varmesh_r[i][j][k+1]<=0))
        {
            if (t_std_varmesh_l[i][j][k] + t_std_varmesh_r[i][j][k+1] > 0)
            {
                t_std_ustar_varmesh[i][j][k]=t_std_varmesh_l[i][j][k];
            }
            else if (t_std_varmesh_l[i][j][k] + t_std_varmesh_r[i][j][k+1] < 0)
            {
                t_std_ustar_varmesh[i][j][k]=t_std_varmesh_r[i][j][k+1];
            }
        }
        else if ((t_std_varmesh_l[i][j][k]<0) && (t_std_varmesh_r[i][j][k+1]>0))
        {
            t_std_ustar_varmesh[i][j][k]=0;
        };
    }
        t_std_ustar_varmesh[i][j][t_std_ustar_varmesh[0].size() -1]=0;
    }


    // Calculate f(u^*)
    // viennacl::matrix<NumericT> t_vie_ustar;
    // viennacl::copy(t_std_ustar_varmesh, t_vie_ustar);
    viennapde::Varmesh<NumericT> t_vie_ustar{t_std_ustar_varmesh};
    // FIXME I have not yet fixed the following part
    // viennacl::matrix<NumericT> t_vie_fustar = viennacl::linalg::element_prod(t_vie_ustar, t_vie_ustar) / 2.0;
    // viennapde::Varmesh<NumericT> t_vie_fustar = t_vie_ustar * t_vie_ustar / 2.0;
    viennapde::Varmesh<NumericT> t_vie_fustar{t_vie_ustar, true};
    t_vie_fustar = t_vie_ustar * t_vie_ustar / 2.0;
    // // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    viennapde::Varmesh<NumericT> t_vie_fustardiff{t_vie_fustar};
    // viennacl::range submat_row_range_from1(0, 1),
    //                 submat_column_range_from1(0, t_vie_fustar.size2()-1),
    //                 submat_row_range_to1(0, 1),
    //                 submat_column_range_to1(1, t_vie_fustardiff.size2());
    // viennacl::matrix_range<viennacl::matrix<NumericT>>  
    //     submatrix_to_minus(t_vie_fustar, submat_row_range_from1, submat_column_range_from1),
    //     submatrix_tobe_minus(t_vie_fustardiff, submat_row_range_to1, submat_column_range_to1);
    // submatrix_tobe_minus -= submatrix_to_minus;
    viennacl::matrix<NumericT> tKernel{3, 3};
    tKernel(1, 0) = -1;
    viennapde::convolve<NumericT, ConvolutionType::EQUIV>(t_vie_fustar, tKernel, t_vie_fustardiff, false);
    // // Time Forward on tVarmesh, attention this is a temporarlly extended Varmesh to suffice the need of periodic B.C..
    tVarmesh -= t_vie_fustardiff * dt / dx;
    // // Go back to the true output Varmesh, we need to cur down the margin.
    // viennacl::range submat_row_range_from2(0, 1),
    //                 submat_column_range_from2(shift, iVarmesh.get_column_num()+shift-1),
    //                 submat_row_range_to2(0, 1),
    //                 submat_column_range_to2(0, iVarmesh.get_column_num());
    // viennacl::matrix_range<viennacl::matrix<NumericT>>  
    //     submatrix_to_back(tVarmesh.data_->at(0), submat_row_range_from2, submat_column_range_from2),
    //     submatrix_tobe_back(oVarmesh.data_->at(0), submat_row_range_to2, submat_column_range_to2);
    // oVarmesh.data_->at(0) = submatrix_to_back;

}

enum Scheme {
    GodunovENUM
};


template <typename NumericT>
void RK3order(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    viennapde::Varmesh<NumericT> & oVarmesh, NumericT dt, NumericT dx, Scheme scheme)
{
    if (scheme==Scheme::GodunovENUM)
    {
        viennapde::Varmesh<NumericT> tVarmesh1(iVarmesh), tVarmesh2(iVarmesh);

        Godunov(iVarmesh, tVarmesh1, dt, dx);
        Godunov(tVarmesh1, tVarmesh2, dt, dx);
        for (size_t i = 0; i < iVarmesh.get_var_num(); i++)
        {
            viennacl::matrix<NumericT> t_matrix(tVarmesh2.data_[i]);
            tVarmesh2.data_[i] = 3 / 4 * iVarmesh.data_[i] + 1 / 4 * t_matrix;
        }
        Godunov(tVarmesh2, oVarmesh, dt, dx);
        for (size_t i = 0; i < iVarmesh.get_var_num(); i++)
        {
            viennacl::matrix<NumericT> t_matrix(oVarmesh.data_[i]);
            oVarmesh.data_[i] = 1 / 3 * iVarmesh.data_[i] + 2 / 3 * t_matrix;
        }
        
    }

}

} //namespace viennapde::scheme


} //namespace viennapde


