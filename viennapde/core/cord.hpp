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

/** @file viennapde/core/convol.hpp
    @brief Implementation of the variable container for PDE numerical solver by matrix class of ViennaCL library.
    The 2D-style data handling utilizes the potential of GPU to its extent.
*/

#include "viennacl/matrix.hpp"




template <typename IntegerT>
struct cord2
{
    IntegerT x;
    IntegerT y;
    explicit cord2(IntegerT iX, IntegerT iY): x{iX}, y{iY} { };
};

template <typename IntegerT>
struct cord3
{
    IntegerT x;
    IntegerT y;
    IntegerT z = 0;
    explicit cord3(IntegerT iX, IntegerT iY, IntegerT iZ): x{iX}, y{iY}, z{iZ} { };
    explicit cord3(IntegerT iX, IntegerT iY): cord3{iX, iY, 0} { };
};

namespace viennapde{
    using GridIntT = ssize_t;
}