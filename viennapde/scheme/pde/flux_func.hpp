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

#include "viennapde/core/mesh.hpp"

namespace viennapde {
namespace scheme {

// TODO  Add more functions to make it generic for various function 
template <typename NumericT, FluxType fluxT>
inline mesh<NumericT> FluxFunction(const mesh<NumericT> & u_neg, const mesh<NumericT> & u_pos, NumericT dt=0, NumericT dx=0)
{
    if constexpr (fluxT == FluxType::LaxFriedrichs_flux) {
        /* LaxFriedrichs Flux function 
        * \hat{f}_{j+1/2} = \hat{f}(u_{j}, u_{j+1}) 
        * = 1/2 [f(u_j) + f(u_{j+1}) - \alpha (u_{j+1}-u_j)], where \alpha = \max_u |f'(u)|, in the Burger's equation, \alpha = \max(u_j,u_{j+1})
        */
        return (u_neg*u_neg + u_pos*u_pos)/4.0 - elem_max(elem_abs(u_neg), elem_abs(u_pos)) * (u_pos-u_neg);
    } else if constexpr (fluxT == FluxType::LaxWendroff_flux) {
        /* LaxWendroff Flux function 
        * \hat{f}_{j+1/2} = \hat{f}(u_{j}, u_{j+1}) 
        * = 1/2 [f(u_j) + f(u_{j+1})] - dt/2dx f'((u_j+u_{j+1})/2)[f(u_{j+1})-f(u_j)] 
        * Here, due to FVM, u_{j} and u_{j+1} are taken by average value \bar{u}_j and \bar{u}_{j+1}.
        */
        return ((u_neg^2) + (u_pos^2))/4.0 - (dt/dx/2 * 0.5 * 0.5) * (u_neg + u_pos) * ((u_pos^2) - (u_neg^2));
    }
}

} //namespace viennapde::scheme
} //namespace viennapde
