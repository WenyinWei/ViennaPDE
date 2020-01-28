#pragma once
/* =========================================================================
   Copyright (c) 2016-2020, Department of Engineering Physics,
                            Tsinghua University, Beijing, China.

   Portions of this software are copyright by UChicago Argonne, LLC and ViennaCL team.

                            -----------------
                  viennapde - The Vienna PDE Library
                            -----------------

   Project Head:    Wenyin Wei                   weiwy16@mails.tsinghua.edu.cn

   (A list of authors and contributors can be found in the manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennapde/core/viannapde_enum.hpp
    @brief Offer constant and enumeration definitions in viennapde library
*/


namespace viennapde
{

enum OptimizeLevel {
    First = 0,
    Second = -1,
    Third = -2,
    Fourth = -3,
};



enum ConvolutionType
{
    INNER,
    OUTER,
    EQUIV
};


} //namespace viennapde


