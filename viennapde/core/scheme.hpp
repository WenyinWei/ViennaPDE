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
#define _USE_MATH_DEFINES
#include <cmath>
#include <functional>

#include "viennacl/linalg/matrix_operations.hpp"
#include "viennapde/core/mesh.hpp"
#include "viennapde/core/meshb.hpp"
#include "viennapde/core/convol_mesh.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
// TODO: Change it to operator overloading
template <typename NumericT>
inline mesh<NumericT> elem_geq0(const mesh<NumericT> & iMesh)
{
    mesh<NumericT> tMesh{iMesh.get_size_num()};
    viennacl::matrix<NumericT> scalar_mat = viennacl::scalar_matrix<NumericT>(iMesh.get_row_num(), iMesh.get_column_num(), 0.5);
    for (size_t i = 0; i < iMesh.get_layer_num(); i++)
        *(tMesh[i]) = M_1_PI * viennacl::linalg::element_atan(*(iMesh[i]) * 1000.0) + scalar_mat;
    return tMesh;
}

namespace scheme
{

template <typename NumericT>
mesh<NumericT> Godunov(
    const mesh<NumericT> & iMesh, // This func requires extend_boundary(0,1,1);
    NumericT dt, NumericT dx)
{
    // iMesh.check_boundary(0,1,0);
    // iMesh.BCPeriodic();

    // Calculate the u_{j+1/2}^-, u_{j+1/2}^+
    
    // std::vector<std::vector<std::vector<NumericT>>> t_std_varmesh_l, 
    //                                                 t_std_varmesh_r,
    //                                                 t_std_ustar_varmesh;
    // viennacl::copy(tMesh, t_std_varmesh_l);
    // viennacl::copy(tMesh, t_std_varmesh_r);
    // viennacl::copy(tMesh, t_std_ustar_varmesh);

    // // Calculate the u^*
    // //TODO I have no idea how to operate this part on GPU.
    // for (size_t i = 0; i < tMesh.get_layer_num(); i++)
    // for (size_t j = 0; j < tMesh.get_row_num(); j++)
    // {
    //     for (size_t k = 0; k < tMesh.get_column_num()-1; k++)
    //     {
    //         if ((t_std_varmesh_l[i][j][k]>=0) && (t_std_varmesh_r[i][j][k+1]>=0)) 
    //         {
    //             t_std_ustar_varmesh[i][j][k]=t_std_varmesh_l[i][j][k];
    //         }
    //         else if ((t_std_varmesh_l[i][j][k]<=0) && (t_std_varmesh_r[i][j][k+1]<=0))
    //         {
    //             t_std_ustar_varmesh[i][j][k]=t_std_varmesh_r[i][j][k+1];
    //         }
    //         else if ((t_std_varmesh_l[i][j][k]>=0) && (t_std_varmesh_r[i][j][k+1]<=0))
    //         {
    //             if (t_std_varmesh_l[i][j][k] + t_std_varmesh_r[i][j][k+1] > 0)
    //             {
    //                 t_std_ustar_varmesh[i][j][k]=t_std_varmesh_l[i][j][k];
    //             }
    //             else if (t_std_varmesh_l[i][j][k] + t_std_varmesh_r[i][j][k+1] < 0)
    //             {
    //                 t_std_ustar_varmesh[i][j][k]=t_std_varmesh_r[i][j][k+1];
    //             }
    //         }
    //         else if ((t_std_varmesh_l[i][j][k]<0) && (t_std_varmesh_r[i][j][k+1]>0))
    //         {
    //             t_std_ustar_varmesh[i][j][k]=0;
    //         };
    //     }
    //     t_std_ustar_varmesh[i][j][t_std_ustar_varmesh[0][0].size() -1]=0;
    // }
    
    

    // Calculate f(u^*)
    mesh<NumericT> fustar{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 2) = 1;
    mesh<NumericT>         ur = viennapde::convolve<NumericT, ConvolutionType::EQUIV>(iMesh, tKernel);
    mesh<NumericT> const & ul = iMesh;
    mesh<NumericT> ustar 
        = (ur + (ul-ur) * elem_geq0(ul+ur)) * ((NumericT)1.0- ((NumericT)1.0-elem_geq0(ul))*elem_geq0(ur) );
    fustar = ustar * ustar / 2 ;
    }

    // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    mesh<NumericT> fustardiff{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 0) = -1; tKernel(0, 1) = 1;
    fustardiff = viennapde::convolve<NumericT, ConvolutionType::EQUIV>(fustar, tKernel);
    }


    // Time Forward, attention this is a temporarlly extended mesh to suffice the need of periodic B.C..
    mesh<NumericT> oMesh = iMesh - fustardiff * (dt / dx);
    return oMesh;
    // // Test Output Part Make it a function after a while
    // static int times = 0;
    // std::vector< std::vector< std::vector<NumericT> > >  stl_varmesh;
    // viennacl::copy(tMesh, stl_varmesh);
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
}


template <typename NumericT>
mesh<NumericT> RK3order(
    const mesh<NumericT> & iMesh,
    NumericT dt, NumericT dx, 
    std::function<mesh<NumericT> (const mesh<NumericT> & , NumericT, NumericT )> scheme)
{
    std::vector< mesh<NumericT> > tMesh;
    for (size_t i = 0; i < 3; i++)
        tMesh.emplace_back(iMesh.get_size_num());
    tMesh[0] = scheme(iMesh,    dt, dx);
    tMesh[1] = scheme(tMesh[0], dt, dx);
    tMesh[1] = iMesh * (3/4) + tMesh[1] * (1/4);
    tMesh[2] = scheme(tMesh[1], dt, dx);
    tMesh[2] = iMesh * (1/3) + tMesh[2] * (2/3);
    return tMesh[2];
}

} //namespace viennapde::scheme


} //namespace viennapde


