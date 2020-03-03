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

/** @file viennapde/scheme/pde/Godunov.hpp
    @brief Godunov's every edge is regarded as a Riemman problem.
*/
#include <vector>

#include "viennapde/core/mesh.hpp"
#include "viennapde/core/convol_mesh.hpp"
#include "viennapde/core/elem_op.hpp"

// SECTION 01b Declare the image class
namespace viennapde
{
namespace scheme
{

// This func requires extend_boundary(0,1,1);
template <typename NumericT>
/**
 * @brief Impose 1D scheme Godunov based on the shock and rarefaction analysis on the Burger's equation. Requires (0, 1, 0) or more boundary.
 * @param  {const mesh<NumericT> &} iMesh : 
 * @param  {NumericT} dt                  : 
 * @param  {NumericT} dx                  : 
 * @return {mesh<NumericT>}               : 
 */
mesh<NumericT> Godunov(const mesh<NumericT> & iMesh, NumericT dt, NumericT dx)
{
    // Calculate f(u^*)
    mesh<NumericT> hatf{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 2) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,2}};
    mesh<NumericT>         ur = convolve(iMesh, tKernel, tROIrc_vec);
    mesh<NumericT> const & ul = iMesh;
    mesh<NumericT> ustar 
        = (ur + (ul-ur) * elem_geq0(ul+ur)) * (1.0- (1.0-elem_geq0(ul))*elem_geq0(ur) );
    hatf = ustar * ustar / 2 ;
    }

    // Calculate f_{j+1/2}(u^*) - f_{j-1/2}(u^*)
    mesh<NumericT> hatfdiff{iMesh.get_size_num()};
    {
    viennacl::matrix<NumericT> tKernel{1, 3}; tKernel(0, 0) = -1; tKernel(0, 1) = 1;
    std::vector<cord2<GridIntT>> tROIrc_vec{{0,0}, {0,1}};
    hatfdiff = convolve(hatf, tKernel, tROIrc_vec);
    }

    // Time Forward
    return iMesh - hatfdiff * (dt / dx);
}

} //namespace viennapde::scheme
} //namespace viennapde
