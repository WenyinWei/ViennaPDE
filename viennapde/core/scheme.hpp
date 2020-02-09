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

#include <functional>

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
    viennapde::Varmesh<NumericT> & oVarmesh, NumericT dt, NumericT dx)
{
    viennapde::BVarmesh<NumericT> tVarmesh{(DequeMat<NumericT>)iVarmesh};
    tVarmesh.extend_boundary(0,1,0);
    tVarmesh.BCPeriodic();

    // Calculate the u_{j+1/2}^-, u_{j+1/2}^+
    std::vector<std::vector<std::vector<NumericT>>> t_std_varmesh_l, 
                                                    t_std_varmesh_r,
                                                    t_std_ustar_varmesh;
    viennacl::copy(tVarmesh, t_std_varmesh_l);
    viennacl::copy(tVarmesh, t_std_varmesh_r);
    viennacl::copy(tVarmesh, t_std_ustar_varmesh);

    // Calculate the u^*
    //TODO I have no idea how to operate this part on GPU.
    for (size_t i = 0; i < tVarmesh.get_layer_num(); i++)
    for (size_t j = 0; j < tVarmesh.get_row_num(); j++)
    {
        for (size_t k = 0; k < tVarmesh.get_column_num()-1; k++)
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
        t_std_ustar_varmesh[i][j][t_std_ustar_varmesh[0][0].size() -1]=0;
    }

    // Calculate f(u^*)
    viennapde::Varmesh<NumericT> t_vie_ustar{t_std_ustar_varmesh};
    viennapde::Varmesh<NumericT> t_vie_fustar{t_vie_ustar.get_size_num()};
    t_vie_fustar = t_vie_ustar * t_vie_ustar / 2 ;

    // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    viennapde::Varmesh<NumericT> t_vie_fustardiff{t_vie_fustar};

    {
    viennacl::matrix<NumericT> tKernel{1, 3};
    tKernel(0, 0) = -1;
    viennapde::convolve<NumericT, ConvolutionType::EQUIV>(t_vie_fustar, tKernel, t_vie_fustardiff, ClrOut::NO);
    }


    // Time Forward on tVarmesh, attention this is a temporarlly extended Varmesh to suffice the need of periodic B.C..
    tVarmesh -= t_vie_fustardiff * (dt / dx);
    // // Test Output Part Make it a function after a while
    // static int times = 0;
    // std::vector< std::vector< std::vector<NumericT> > >  stl_varmesh;
    // viennacl::copy(tVarmesh, stl_varmesh);
    // std::ofstream file;
    // file.open(std::to_string(times++) + "_tVarmesh.csv");
    // if (file.is_open())
    // {
    //   file << stl_varmesh[0][0][0];        
    //   for (size_t column_i = 1; column_i < stl_varmesh[0][0].size(); column_i++)
    //   {
    //     file << ", "<< stl_varmesh[0][0][column_i];
    //   }
    //   file.close();
    // }

    
    { // STUB Cut down the margin.
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 1) = 1;
    viennapde::convolve<NumericT, ConvolutionType::INNER>(tVarmesh, tKernel, oVarmesh, ClrOut::YES);
    }
}


template <typename NumericT>
void RK3order(
    const viennapde::Varmesh<NumericT> & iVarmesh,
    viennapde::Varmesh<NumericT> & oVarmesh, NumericT dt, NumericT dx, 
    std::function<void (const Varmesh<NumericT> & , Varmesh<NumericT> & , NumericT, NumericT )> scheme)
{
    std::vector<viennapde::Varmesh<NumericT>> tVarmesh;
    for (size_t i = 0; i < 2; i++)
        tVarmesh.push_back(viennapde::Varmesh<NumericT>{iVarmesh.get_size_num()});
    scheme(iVarmesh, tVarmesh[0], dt, dx);
    scheme(tVarmesh[0], tVarmesh[1], dt, dx);
    tVarmesh[1] = iVarmesh * (3/4) + tVarmesh[1] * (1/4);
    scheme(tVarmesh[1], oVarmesh, dt, dx);
    oVarmesh = iVarmesh * (1/3) + oVarmesh * (2/3);
}

} //namespace viennapde::scheme


} //namespace viennapde


