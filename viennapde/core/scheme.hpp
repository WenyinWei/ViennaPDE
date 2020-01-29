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
#include "viennapde/core/varmesh.hpp"


// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

template <typename NumericT>
void Godunov(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    viennapde::Varmesh<NumericT> & o_varmesh, NumericT dt, NumericT dx, unsigned int order=1)
{
    ssize_t shift = 1;
    viennapde::Varmesh<NumericT> t_varmesh(1, iVarmesh.get_row_num(), iVarmesh.get_column_num()+2*shift);
    viennacl::range submat_row_range_from(0, 1),
                    submat_column_range_from(0, iVarmesh.get_column_num()),
                    submat_row_range_to(0, 1),
                    submat_column_range_to(shift, iVarmesh.get_column_num()+shift);
    viennacl::matrix_range<viennacl::matrix<NumericT>>  
        submatrix_to_replace(iVarmesh.data_->at(0), submat_row_range_from, submat_column_range_from),
        submatrix_tobe_replace(t_varmesh.data_->at(0), submat_row_range_to, submat_column_range_to); 
    submatrix_tobe_replace = submatrix_to_replace;
    for (ssize_t i = 0; i < shift; i++)
    {
        t_varmesh.data_->at(0)(0, i) = iVarmesh.data_->at(0)(0, iVarmesh.get_column_num()+i);
        t_varmesh.data_->at(0)(0, iVarmesh.get_column_num()+2*shift-1-i) = iVarmesh.data_->at(0)(0, 2*shift-1-i);
    }
    


    // Calculate the u_{j+1/2}^-, u_{j+1/2}^+
    std::vector<std::vector<NumericT>> t_std_varmesh_l(1), t_std_varmesh_r(1);
    std::vector<std::vector<NumericT>> t_std_ustar_varmesh(1);
    t_std_varmesh_l[0].resize(t_varmesh.get_column_num());
    t_std_varmesh_r[0].resize(t_varmesh.get_column_num());
    t_std_ustar_varmesh[0].resize(t_varmesh.get_column_num());
    
    viennacl::copy(t_varmesh.data_->at(0), t_std_varmesh_l);
    viennacl::copy(t_varmesh.data_->at(0), t_std_varmesh_r);
    viennacl::copy(t_varmesh.data_->at(0), t_std_ustar_varmesh);


    // Calculate the u^*
    for (ssize_t i = 0; i < t_varmesh.get_column_num()-1; i++)
    {
        if ((t_std_varmesh_l[0][i]>=0) && (t_std_varmesh_r[0][i+1]>=0)) 
        {
            t_std_ustar_varmesh[0][i]=t_std_varmesh_l[0][i];
        }
        else if ((t_std_varmesh_l[0][i]<=0) && (t_std_varmesh_r[0][i+1]<=0))
        {
            t_std_ustar_varmesh[0][i]=t_std_varmesh_r[0][i+1];
        }
        else if ((t_std_varmesh_l[0][i]>=0) && (t_std_varmesh_r[0][i+1]<=0))
        {
            if (t_std_varmesh_l[0][i] + t_std_varmesh_r[0][i+1] > 0)
            {
                t_std_ustar_varmesh[0][i]=t_std_varmesh_l[0][i];
            }
            else if (t_std_varmesh_l[0][i] + t_std_varmesh_r[0][i+1] < 0)
            {
                t_std_ustar_varmesh[0][i]=t_std_varmesh_r[0][i+1];
            }
        }
        else if ((t_std_varmesh_l[0][i]<0) && (t_std_varmesh_r[0][i+1]>0))
        {
            t_std_ustar_varmesh[0][i]=0;
        };
    }
    t_std_ustar_varmesh[0][t_std_ustar_varmesh[0].size() -1]=0;

    // Calculate f(u^*)
    viennacl::matrix<NumericT> t_vie_ustar;
    viennacl::copy(t_std_ustar_varmesh, t_vie_ustar);
    viennacl::matrix<NumericT> t_vie_fustar = viennacl::linalg::element_prod(t_vie_ustar, t_vie_ustar) / 2.0;
    // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    viennacl::matrix<NumericT> t_vie_fustardiff = t_vie_fustar;
    viennacl::range submat_row_range_from1(0, 1),
                    submat_column_range_from1(0, t_vie_fustar.size2()-1),
                    submat_row_range_to1(0, 1),
                    submat_column_range_to1(1, t_vie_fustardiff.size2());
    viennacl::matrix_range<viennacl::matrix<NumericT>>  
        submatrix_to_minus(t_vie_fustar, submat_row_range_from1, submat_column_range_from1),
        submatrix_tobe_minus(t_vie_fustardiff, submat_row_range_to1, submat_column_range_to1);
    submatrix_tobe_minus -= submatrix_to_minus;
    // Time Forward on t_varmesh, attention this is a temporarlly extended Varmesh to suffice the need of periodic B.C..
    t_varmesh.data_->at(0) -= dt / dx * t_vie_fustardiff;

    // Go back to the true output Varmesh, we need to cur down the margin.
    viennacl::range submat_row_range_from2(0, 1),
                    submat_column_range_from2(shift, iVarmesh.get_column_num()+shift-1),
                    submat_row_range_to2(0, 1),
                    submat_column_range_to2(0, iVarmesh.get_column_num());
    viennacl::matrix_range<viennacl::matrix<NumericT>>  
        submatrix_to_back(t_varmesh.data_->at(0), submat_row_range_from2, submat_column_range_from2),
        submatrix_tobe_back(o_varmesh.data_->at(0), submat_row_range_to2, submat_column_range_to2);
    o_varmesh.data_->at(0) = submatrix_to_back;

}

enum Scheme {
    GodunovENUM
};


template <typename NumericT>
void RK3order(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    viennapde::Varmesh<NumericT> & o_varmesh, NumericT dt, NumericT dx, Scheme scheme)
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
        Godunov(tVarmesh2, o_varmesh, dt, dx);
        for (size_t i = 0; i < iVarmesh.get_var_num(); i++)
        {
            viennacl::matrix<NumericT> t_matrix(o_varmesh.data_[i]);
            o_varmesh.data_[i] = 1 / 3 * iVarmesh.data_[i] + 2 / 3 * t_matrix;
        }
        
    }

}

} //namespace viennapde::scheme


} //namespace viennapde


