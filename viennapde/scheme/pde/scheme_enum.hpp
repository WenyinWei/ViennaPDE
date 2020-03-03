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

/** @file viennapde/scheme/scheme_enum.hpp
    @brief Implementation of concrete scheme for PDE.
*/

// SECTION 01b Declare the image class
namespace viennapde {
namespace scheme {

enum FluxType {
    LaxFriedrichs_flux,
    LaxWendroff_flux,
    Godunov_flux
};

} //namespace viennapde::scheme
} //namespace viennapde
