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
#ifndef _HEARTDIFFUSIONFUNCTOR_H_
#define _HEARTDIFFUSIONFUNCTOR_H_
#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <lifev/core/filter/GetPot.hpp>
#include <lifev/core/fem/FESpace.hpp>
#include "StimElectrode.hpp"

namespace LifeV
{

class HeartDiffusionFunctor
{
public:
    typedef boost::function<Real ( Real const& x, Real const& y, Real const& z, Real const&, ID const& id) > funcType;
    typedef std::shared_ptr<ripple::StimElectrode<Real, Real>> elecPtrType;
    typedef ripple::AxisStimElectrode<Real, Real> axisElecType;
    typedef std::shared_ptr<axisElecType> axisElecPtrType;
private:
    GetPot& M_dataFile;
    std::string M_stimulusMode;
    Real M_restPotential;
    std::vector<elecPtrType> vecPtrStimElectrode;
private:
    HeartDiffusionFunctor ( const HeartDiffusionFunctor& heartFunctors );
    HeartDiffusionFunctor& operator= ( const HeartDiffusionFunctor& heartFunctors );
    Real setStimulus ( const Real& t, const Real& x, const Real& y, const Real& z, const ID&   id) const;
    Real setInitialScalar ( const Real& t, const Real& x, const Real& y, const Real& z, const ID& id );
    Real setZeroScalar ( const Real& t, const Real& x, const Real& y, const Real& z , const ID& id );
public:
    HeartDiffusionFunctor();
    HeartDiffusionFunctor ( GetPot& dataFile );
    virtual ~HeartDiffusionFunctor(){};
    funcType stimulus();
    funcType initialScalar();
    funcType zeroScalar();
};

} // namespace LifeV
#endif
