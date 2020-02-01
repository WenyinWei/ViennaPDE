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

/** @file viennapde/core/mvarmesh.hpp
    @brief Varmesh with margin.
*/

#include "./varmesh.hpp"




// SECTION 01 Define the variable mesh (Varmesh) class
namespace viennapde
{

template <typename NumericT>
class MVarmesh : public Varmesh<NumericT>
{
private: 
    cord3<GridIntT> margin_;
public:
    cord3<GridIntT> get_margin() const { return margin_;};
protected:
    void set_margin(GridIntT iX, GridIntT iY, GridIntT iZ) {margin_.x=iX; margin_.y=iY; margin_.z=iZ;};
    
    
};







} //namespace viennapde


