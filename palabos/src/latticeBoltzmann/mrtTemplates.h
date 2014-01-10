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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef MRT_TEMPLATES_H
#define MRT_TEMPLATES_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/mrtLattices.h"

namespace plb {
    
template<typename T, class Descriptor> struct mrtTemplatesImpl;

/// All helper functions are inside this structure
template<typename T, template<typename U> class Descriptor>
struct mrtTemplates {

    /// Computation of equilibrium distribution (in moments space)
    static T equilibrium( plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                          const T jSqr ) {
        mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>
            ::equilibrium(iPop, rhoBar, j, jSqr );
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor<T>::q>& momentsEq, 
                                    T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                    const T jSqr ) {
        
        mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            computeEquilibrium( momentsEq, rhoBar, j, jSqr );
    }
    
    static void computeMoments(Array<T,Descriptor<T>::q> &moments, Cell<T,Descriptor> &cell)
    {
        mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            computeMoments(moments, cell.getRawPopulations());
    }
    
    /// MRT collision step
    static T mrtCollision( Cell<T,Descriptor>& cell, T omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
        mrtCollision( cell.getRawPopulations(), omega);
    }
    
    /// MRT collision step (imposed rhoBar, j)
    static T mrtCollision( Cell<T,Descriptor>& cell,
                           T const& rhoBar, Array<T,Descriptor<T>::d> const& j,
                           const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            mrtCollision( cell.getRawPopulations(), rhoBar, j, omega);
    }
    
    /// quasi incompressible MRT collision step with force
    static T smagorinskyMrtCollision( Cell<T,Descriptor>& cell,
                                        const T &rhoBar, const Array<T,Descriptor<T>::d> & j,
                                        const Array<T,SymmetricTensor<T,Descriptor>::n > &strain, T cSmago, 
                                        const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor >::
            smagorinskyMrtCollision( cell.getRawPopulations(), rhoBar, j, omega, strain, cSmago);
    }
    
    /// MRT collision step with force
    static T mrtCollisionWithForce( Cell<T,Descriptor>& cell,
                           const T &rhoBar, const Array<T,Descriptor<T>::d> & u,
                           const T &omega, T amplitude)
    {
        Array<T,Descriptor<T>::d> force;
        force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor >::
            mrtCollisionWithForce( cell.getRawPopulations(), rhoBar, u, omega, force,amplitude);
    }
    
    /// MRT collision step with force
    static T smagorinskyMrtCollisionWithForce( Cell<T,Descriptor>& cell,
                                               const T &rhoBar, const Array<T,Descriptor<T>::d> & u,
                                               const Array<T,SymmetricTensor<T,Descriptor>::n > &strain, T cSmago, 
                                               const T &omega, T amplitude)
    {
        Array<T,Descriptor<T>::d> force;
        force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor >::
        smagorinskyMrtCollisionWithForce( cell.getRawPopulations(), rhoBar, u, omega, strain, cSmago, force, amplitude);
    }
    
    
    /// quasi incompressible MRT collision step with force
    static T incMrtCollisionWithForce( Cell<T,Descriptor>& cell,
                                    const T &rhoBar, const Array<T,Descriptor<T>::d> & u,
                                    const T &omega, T amplitude)
    {
        Array<T,Descriptor<T>::d> force;
        force.from_cArray(cell.getExternal(Descriptor<T>::ExternalField::forceBeginsAt));
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor >::
            incMrtCollisionWithForce( cell.getRawPopulations(), rhoBar, u, omega, force,amplitude);
    }
    
    /// quasi incompressible MRT collision step with force
    static T incSmagorinskyMrtCollisionWithForce( Cell<T,Descriptor>& cell,
                                            const T &rhoBar, const Array<T,Descriptor<T>::d> & u,
                                            const Array<T,SymmetricTensor<T,Descriptor>::n > &strain, T cSmago, 
                                            const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor >::
            incSmagorinskyMrtCollisionWithForce( cell.getRawPopulations(), rhoBar, u, omega, strain, cSmago);
    }
    
    /// Add a force term after BGK collision, according to the Guo algorithm
    static void addGuoForce( Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u,
                             const T &omega, T amplitude)
    {
        mrtTemplatesImpl<T, typename Descriptor<T>::SecondBaseDescriptor >
            ::addGuoForce(cell.getRawPopulations(), cell.getExternal(0), u, omega, amplitude);
    }
    
     /// Computation of all equilibrium distribution (in moments space)
    static void computeIncEquilibrium( Array<T,Descriptor<T>::q>& momentsEq, 
                                    T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                    const T jSqr ) {
        
        mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            computeIncEquilibrium( momentsEq, rhoBar, j, jSqr );
    }
    
    /// MRT collision step
    static T incMrtCollision( Cell<T,Descriptor>& cell, const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            incMrtCollision( cell.getRawPopulations(), omega);
    }
    
    /// MRT collision step
    static T incMrtCollision( Cell<T,Descriptor>& cell,
                                   const T &rhoBar, const Array<T,Descriptor<T>::d> & j,
                                   const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            incMrtCollision( cell.getRawPopulations(), rhoBar, j, omega);
    }
    
    /// Computation of all equilibrium distribution (in moments space) for the smagorinsky model
    static void computeIncSmagorinskyEquilibrium( Array<T,Descriptor<T>::q>& momentsEq, 
                                    T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                    const T jSqr,
                                    const Array<T,SymmetricTensor<T,Descriptor>::n > &strain, T cSmago ) 
    {
        
        mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            computeIncSmagorinskyEquilibrium( momentsEq, rhoBar, j, jSqr, strain, cSmago );
    }
    
    /// MRT collision step
    static T incSmagorinskyMrtCollision( Cell<T,Descriptor>& cell,
                                              T &rhoBar, const Array<T,Descriptor<T>::d> & j,
                                              const Array<T,SymmetricTensor<T,Descriptor>::n > &strain, T cSmago,
                                              const T &omega)
    {
        return mrtTemplatesImpl<T,typename Descriptor<T>::SecondBaseDescriptor>::
            incSmagorinskyMrtCollision( cell.getRawPopulations(), rhoBar, j, omega, strain, cSmago);
    }


};  // struct mrtHelpers

template<typename T, class Descriptor>
struct mrtTemplatesImpl {
    
    /// Computation of equilibrium distribution (in moments space)
    static T equilibrium( plint iPop, T rhoBar, Array<T,Descriptor::d> const& j, const T jSqr ) {
        T invRho = Descriptor::invRho(rhoBar);
        T equ = T();
        for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
            equ += Descriptor::M[iPop][jPop] * dynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma2_equilibrium(jPop, rhoBar, invRho, j, jSqr);
        }
        
        return equ;
    }
    
    static void computeInvM_S(T invM_S[Descriptor::q][Descriptor::q], const T &omega) {
        Array<T,Descriptor::q> s;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            s[iPop] = Descriptor::S[iPop];
        }
        for (plint iA = 0; iA < Descriptor::shearIndexes; ++iA) {
            plint iPop = Descriptor::shearViscIndexes[iA];
            s[iPop] = omega;
        }
        
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                invM_S[iPop][jPop] = Descriptor::invM[iPop][jPop] * s[jPop];
            }
        }
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibriumMoments( Array<T,Descriptor::q>& momentsEq, 
                                           T rhoBar, Array<T,Descriptor::d> const& j,
                                           const T jSqr ) {
        T invRho = Descriptor::invRho(rhoBar);
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            momentsEq[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                momentsEq[iPop] += Descriptor::M[iPop][jPop] * 
                    dynamicsTemplatesImpl<T,Descriptor>::
                        bgk_ma2_equilibrium(jPop, rhoBar, invRho, j, jSqr);
            }
        }
    }
    
    static void computeMoments(Array<T,Descriptor::q> &moments, 
                               const Array<T,Descriptor::q>& f) {
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            moments[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                moments[iPop] += Descriptor::M[iPop][jPop] * f[jPop];
            }
        }
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,Descriptor::d> & j,
                           const T &omega) {
        Array<T,Descriptor::q> momentsEq;
        Array<T,Descriptor::q> moments;
        
        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
        computeEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        
        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);
    
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }
        
        return jSqr;
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f,const T &omega) {
        Array<T,Descriptor::q> momentsEq;
        Array<T,Descriptor::q> moments;
        
        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,Descriptor::d> j;
        for (plint iA = 0; iA < Descriptor::jIndexes; ++iA) {
            plint iPop = Descriptor::momentumIndexes[iA];
            
            j[iA] = moments[iPop];
        }
        T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
        computeEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        
        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);
    
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }
        
        return jSqr;
    }
    
    /// smagorinsky MRT collision step
    static T smagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                                      const T &rhoBar, const Array<T,Descriptor::d> & j,
                                      const T &omega, 
                                      const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n >& strain, T cSmago ) 
    {
        PLB_ASSERT(false);
    }
    
    static void addGuoForce( Array<T,Descriptor::q>& f, const Array<T,Descriptor::d> &force,
                             Array<T,Descriptor::d> const& u,
                             const T &omega, T amplitude )
    {
        Array<T,Descriptor::q> forcing;
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T c_u = T();
            for (int iD=0; iD<Descriptor::d; ++iD) {
                c_u += Descriptor::c[iPop][iD]*u[iD];
            }
            c_u *= Descriptor::invCs2 * Descriptor::invCs2;
            T forceTerm = T();
            for (int iD=0; iD < Descriptor::d; ++iD) {
                forceTerm +=
                    (   ((T)Descriptor::c[iPop][iD]-u[iD]) * Descriptor::invCs2
                    + c_u * (T)Descriptor::c[iPop][iD]
                    )
                    * force[iD];
            }
            forceTerm *= Descriptor::t[iPop];
            forceTerm *= amplitude;
            forcing[iPop] = forceTerm;
        }
        
        Array<T,Descriptor::q> forceMoments;
        computeMoments(forceMoments,forcing);
        
        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);
        
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += -invM_S[iPop][jPop] * forceMoments[jPop];
            }
            collisionTerm *= (T)0.5;
            collisionTerm += forcing[iPop];
            f[iPop] += collisionTerm;
        }
    }
    
    /// MRT collision step with force
    static T mrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    const T &omega, 
                                    const Array<T,Descriptor::d>& force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = mrtCollision( f, rhoBar, j, omega );
        addGuoForce( f, force, u, omega, amplitude );
        
        return jSqr;
    }
    
    /// MRT collision step with force
    static T smagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                               const T &rhoBar, const Array<T,Descriptor::d> & u,
                                               const T &omega, 
                                               const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n >& strain, T cSmago,
                                               const Array<T,Descriptor::d>& force, T amplitude) 
    {
        PLB_ASSERT(false);
    }
    
    /// quasi incompressible MRT collision step with force
    static T incMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    const T &omega, 
                                    const Array<T,Descriptor::d>& force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = incMrtCollision( f, rhoBar, j, omega);
        addGuoForce( f, force, u, omega, amplitude );
        
        return jSqr;
    }
    
    /// quasi incompressible MRT collision step with force
    static T incSmagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                                       const T &rhoBar, const Array<T,Descriptor::d> & u,
                                                       const T &omega, 
                                                       const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n >& strain, T cSmago,
                                                       const Array<T,Descriptor::d>& force, T amplitude ) 
    {
        PLB_ASSERT(false);
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncEquilibrium( Array<T,Descriptor::q>& momentsEq, 
                                            T rhoBar, Array<T,Descriptor::d> const& j,
                                            const T jSqr )
    {
        T invRho = (T)1;
        for (plint iPop = 0; iPop < Descriptor::q; ++iPop) {
            momentsEq[iPop] = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                momentsEq[iPop] += Descriptor::M[iPop][jPop] * 
                    dynamicsTemplatesImpl<T,Descriptor>::
                        bgk_ma2_equilibrium(jPop, rhoBar, invRho, j, jSqr);
            }
        }
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncSmagorinskyEquilibrium( Array<T,Descriptor::q>& momentsEq, 
                                            T rhoBar, Array<T,Descriptor::d> const& j,
                                            const T jSqr,
                                            const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n >& strain, T cSmago )
    {
        PLB_ASSERT(false);
    }
    
    /// MRT collision step
    static T incMrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,Descriptor::d> & j,
                           const T &omega) {
        Array<T,Descriptor::q> momentsEq;
        Array<T,Descriptor::q> moments;
        
        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
        computeIncEquilibriumMOments(momentsEq,rhoBar,j,jSqr);
        
        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);
    
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }
        
        return jSqr;
    }
    
    /// MRT collision step
    static T incMrtCollision( Array<T,Descriptor::q>& f,
                           const T &omega) {
        Array<T,Descriptor::q> momentsEq;
        Array<T,Descriptor::q> moments;
        
        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,Descriptor::d> j;
        for (plint iA = 0; iA < Descriptor::jIndexes; ++iA) {
            plint iPop = Descriptor::momentumIndexes[iA];
            
            j[iA] = moments[iPop];
        }
        T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
        computeIncEquilibriumMOments(momentsEq,rhoBar,j,jSqr);
        
        T invM_S[Descriptor::q][Descriptor::q];
        computeInvM_S(invM_S, omega);
    
        for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
            T collisionTerm = T();
            for (plint jPop = 0; jPop < Descriptor::q; ++jPop) {
                collisionTerm += invM_S[iPop][jPop] * (moments[jPop] - momentsEq[jPop]);
            }
            f[iPop] -= collisionTerm;
        }
        
        return jSqr;
    }
    
    /// Smagorinsky MRT collision step
    static T incSmagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                                   const T &rhoBar, const Array<T,Descriptor::d> & j,
                                   const T &omega,
                                   const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n >& strain, T cSmago) 
    {
        PLB_ASSERT(false);
    }
};

}  // namespace plb

// The specialized code is directly included. That is because we never want
// it to be precompiled so that in both the precompiled and the
// "include-everything" version, the compiler can apply all the
// optimizations it wants.
#include "latticeBoltzmann/mrtTemplates2D.h"
#include "latticeBoltzmann/mrtTemplates3D.h"

#endif
