/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OFF_LATTICE_BOUNDARY_PROFILES_3D_H
#define OFF_LATTICE_BOUNDARY_PROFILES_3D_H

#include "core/globalDefs.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "algorithm/functions.h"

namespace plb {


template<typename T, class SurfaceData>
class BoundaryProfile3D {
public:
    virtual ~BoundaryProfile3D() { }
    virtual void setNormal(Array<T,3> const& normal_) =0;
    virtual void defineCircularShape(Array<T,3> const& radius_, T center_) =0;
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          SurfaceData& data, OffBoundary::Type& bdType ) const =0;
    virtual BoundaryProfile3D<T,SurfaceData>* clone() const =0;
};

template<typename T, class SurfaceData>
struct DefaultWallProfile3D {
    BoundaryProfile3D<T,SurfaceData>* generate() {
        // A default wall profile needs yet to be implemented for this case.
        PLB_ASSERT( false );
        return 0;
    }
}; 

template<typename T, class SurfaceData>
BoundaryProfile3D<T, SurfaceData>* generateDefaultWallProfile3D() {
    return DefaultWallProfile3D<T,SurfaceData>().generate();
}

template<typename T>
class NoSlipProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual NoSlipProfile3D<T>* clone() const;
};

template<typename T>
class FreeSlipProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual FreeSlipProfile3D<T>* clone() const;
};

template<typename T>
class ConstantVelocityProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    ConstantVelocityProfile3D(Array<T,3> const& u_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual ConstantVelocityProfile3D<T>* clone() const;
private:
    Array<T,3> u;
};

template<typename T>
class VelocityPlugProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    VelocityPlugProfile3D(T uMax_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual VelocityPlugProfile3D<T>* clone() const;
private:
    T uMax;
    Array<T,3> normal;
};

template< typename T, template<typename U> class Descriptor>
class OscillatingPoiseuilleProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    OscillatingPoiseuilleProfile3D(T minUave_, T maxUave_, T period_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual OscillatingPoiseuilleProfile3D<T,Descriptor>* clone() const;
private:
    T minUave, maxUave, period;
    Array<T,3> normal;
    Array<T,3> center;
    T radius;
};

template< typename T, template<typename U> class Descriptor>
class IncreasingPoiseuilleProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    IncreasingPoiseuilleProfile3D(T uAverage_, plint maxT_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual IncreasingPoiseuilleProfile3D<T,Descriptor>* clone() const;
private:
    T uAverage;
    plint maxT;
    Array<T,3> normal;
    Array<T,3> center;
    T radius;
};

template< typename T, template<typename U> class Descriptor>
class IncreasingVelocityProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    IncreasingVelocityProfile3D(Array<T,3> const& u_, plint maxT_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual IncreasingVelocityProfile3D<T,Descriptor>* clone() const;
private:
    Array<T,3> u;
    plint maxT;
};

template< typename T, template<typename U> class Descriptor>
class TimeDependentVelocityProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    TimeDependentVelocityProfile3D(util::TimeDependentFunction<T,3>* velocity_);
    TimeDependentVelocityProfile3D(TimeDependentVelocityProfile3D<T,Descriptor> const& rhs);
    TimeDependentVelocityProfile3D<T,Descriptor>& operator= (
            TimeDependentVelocityProfile3D<T,Descriptor> const& rhs );
    void swap(TimeDependentVelocityProfile3D<T,Descriptor>& rhs);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual TimeDependentVelocityProfile3D<T,Descriptor>* clone() const;
private:
    util::TimeDependentFunction<T,3>* velocity;
};

template<typename T>
class PoiseuilleProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    PoiseuilleProfile3D(T uAverage_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual PoiseuilleProfile3D<T>* clone() const;
private:
    T uAverage;
    Array<T,3> normal;
    Array<T,3> center;
    T radius;
};

template<typename T>
class NeumannBoundaryProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual NeumannBoundaryProfile3D<T>* clone() const;
private:
    Array<T,3> normal;
};

template<typename T>
class DensityNeumannBoundaryProfile3D : public BoundaryProfile3D<T, Array<T,3> >
{
public:
    DensityNeumannBoundaryProfile3D(T rho_ = (T)1);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,3>& data, OffBoundary::Type& bdType ) const;
    virtual DensityNeumannBoundaryProfile3D<T>* clone() const;
private:
    Array<T,3> normal;
    T rho;
};

template<typename T>
class ScalarNeumannProfile3D : public BoundaryProfile3D<T,Array<T,2> >
{
public:
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarNeumannProfile3D<T>* clone() const;
private:
    Array<T,3> normal;
};

template<typename T>
class ScalarDirichletProfile3D : public BoundaryProfile3D<T,Array<T,2> >
{
public:
    ScalarDirichletProfile3D(T value_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarDirichletProfile3D<T>* clone() const;
private:
    T value;
};

template<typename T>
class ScalarFluxProfile3D : public BoundaryProfile3D<T,Array<T,2> >
{
public:
    ScalarFluxProfile3D(T gradVal_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarFluxProfile3D<T>* clone() const;
private:
    Array<T,3> normal;
    T gradVal;
};

/// Implements the condition grad(rho) = kappa(asymptoticRho-rho).
template<typename T>
class ScalarIsolationProfile3D : public BoundaryProfile3D<T,Array<T,2> >
{
public:
    ScalarIsolationProfile3D(T asymptoticRho_, T kappa_);
    virtual void setNormal(Array<T,3> const& normal_);
    virtual void defineCircularShape(Array<T,3> const& center_, T radius_);
    virtual void getData( Array<T,3> const& pos, plint id, AtomicBlock3D const* argument,
                          Array<T,2>& data, OffBoundary::Type& bdType ) const;
    virtual ScalarIsolationProfile3D<T>* clone() const;
private:
    Array<T,3> normal;
    T asymptoticRho, kappa;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_PROFILES_3D_H

