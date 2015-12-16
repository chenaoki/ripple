//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

#include "HeartDiffusionFunctor.hpp"
#include <boost/shared_ptr.hpp>

using namespace LifeV;
// ===================================================
// Set Methods
// ===================================================

HeartDiffusionFunctor::HeartDiffusionFunctor(  ):
    M_dataFile(*(boost::shared_ptr<GetPot>( new GetPot() ) ) ),
    M_stimulusMode(""),
    M_restPotential(-84.)
{}

HeartDiffusionFunctor::HeartDiffusionFunctor ( GetPot& dataFile ):
    M_dataFile(dataFile),
    M_stimulusMode("S1S2"),
    M_restPotential(dataFile ("electric/physics/u0", -84.))
{
    M_stimulusMode = M_dataFile("", "S1S2");
}

Real
HeartDiffusionFunctor::setStimulus ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const
{
    Real ret = 0.;
    if(M_stimulusMode.compare("file")){}
    return ret;
}

Real
HeartDiffusionFunctor::setInitialScalar (
    const Real& /*t*/,
    const Real& /*x*/,
    const Real& /*y*/,
    const Real& /*z*/,
    const ID&   /*id*/ )
{
    return M_restPotential;
}



Real HeartDiffusionFunctor::setZeroScalar (
    const Real& /*t*/,
    const Real& /*x*/,
    const Real& /*y*/,
    const Real& /*z*/,
    const ID&   /*id*/ )
{
    return 0.;
}

// ===================================================
// Get Methods
// ===================================================

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::stimulus()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setStimulus, this, _1, _2, _3, _4, _5);
    return f;
}

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::initialScalar()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setInitialScalar, this, _1, _2, _3, _4, _5);
     return f;
}

HeartDiffusionFunctor::funcType
HeartDiffusionFunctor::zeroScalar()
{
    funcType f;
    f =  boost::bind (&HeartDiffusionFunctor::setZeroScalar, this, _1, _2, _3, _4, _5);
    return f;
}
