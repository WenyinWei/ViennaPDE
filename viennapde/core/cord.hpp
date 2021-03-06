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

/** @file viennapde/core/cord.hpp
    @brief Structures of cordinate to facilitate cord passing.
*/



template <typename IntegerT>
struct cord2
{
    IntegerT x;
    IntegerT y;
    cord2(IntegerT iX, IntegerT iY): x{iX}, y{iY} { };
};

template <typename IntegerT>
struct cord3
{
    IntegerT x;
    IntegerT y;
    IntegerT z = 0;
    cord3(IntegerT iX, IntegerT iY, IntegerT iZ): x{iX}, y{iY}, z{iZ} { };
    cord3(IntegerT iX, IntegerT iY): cord3{iX, iY, 0} { };
};

template <typename IntegerT>
inline bool operator== (const cord2<IntegerT> & lhs, const cord2<IntegerT> & rhs) { return (lhs.x==rhs.x && lhs.y==rhs.y);}
template <typename IntegerT>
inline bool operator== (const cord3<IntegerT> & lhs, const cord3<IntegerT> & rhs) { return (lhs.x==rhs.x && lhs.y==rhs.y && lhs.z==rhs.z);}