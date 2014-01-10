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
 * Template specializations for some computationally intensive LB 
 * functions of the header file mrtTemplates.h, for the D3Q19 grid.
 */

#ifndef MRT_TEMPLATES_3D_H
#define MRT_TEMPLATES_3D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D3Q19 lattice
template<typename T> struct mrtTemplatesImpl<T, descriptors::MRTD3Q19DescriptorBase<T> > {
    
    typedef descriptors::D3Q19DescriptorBase<T> Descriptor; 
    typedef descriptors::MRTD3Q19DescriptorBase<T> MRTDescriptor; 

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibriumMoments( Array<T,Descriptor::q>& momentsEq,
                                    T rhoBar, Array<T,3> const& j, T jSqr )
    {
        T invRho = Descriptor::invRho(rhoBar);
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*jSqr*invRho-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*jSqr*invRho+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*j[0]*j[0]-j[1]*j[1]-j[2]*j[2])*invRho;
        momentsEq[10] = (-j[0]*j[0]+(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2])*invRho;
        momentsEq[11] = (j[1]*j[1]-j[2]*j[2])*invRho;
        momentsEq[12] = (-(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2])*invRho;
        momentsEq[13] = j[1]*j[0]*invRho;
        momentsEq[14] = j[2]*j[1]*invRho;
        momentsEq[15] = j[2]*j[0]*invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
        
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeSmagorinskyEquilibriumMoments( Array<T,Descriptor::q>& momentsEq,
                                               T rhoBar, Array<T,3> const& j, T jSqr, const Array<T,6> &strain, T cSmago )
    {
        typedef SymmetricTensorImpl<T,3> S;
        
        T invRho = Descriptor::invRho(rhoBar);
        T rho = Descriptor::fullRho(rhoBar);
        T rho2 = rho*rho;
        
        T sNorm = sqrt((T)2*SymmetricTensorImpl<T,3>::tensorNormSqr(strain));
        T smagoFactor = (T)2*cSmago * cSmago * sNorm;
        
        T ux2 = rho2 * smagoFactor * strain[S::xx];
        T uy2 = rho2 * smagoFactor * strain[S::yy];
        T uz2 = rho2 * smagoFactor * strain[S::zz];
        
        T uxuy = rho2 * smagoFactor * strain[S::xy];
        T uyuz = rho2 * smagoFactor * strain[S::yz];
        T uxuz = rho2 * smagoFactor * strain[S::xz];
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*(jSqr + ux2+uy2+uz2)*invRho-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*(jSqr + ux2+uy2+uz2)*invRho + (T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*(j[0]*j[0]+ux2)-(j[1]*j[1]+uy2)-(j[2]*j[2]+uz2))*invRho;
        momentsEq[10] = (-(j[0]*j[0]+ux2)+(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2))*invRho;
        momentsEq[11] = (j[1]*j[1]+uy2-(j[2]*j[2]+uz2))*invRho;
        momentsEq[12] = (-(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2))*invRho;
        momentsEq[13] = (j[1]*j[0]+uxuy)*invRho;
        momentsEq[14] = (j[2]*j[1]+uyuz)*invRho;
        momentsEq[15] = (j[2]*j[0]+uxuz)*invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// Computation of all moments (specialized for d3q19)
    static void computeMoments(Array<T,Descriptor::q>& moments, const Array<T,Descriptor::q>& f)
    {
        moments[0] = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[10]+f[11]+f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
        moments[1] =
            -(T)30*f[0]-(T)11*(f[1]+f[2]+f[3]+f[10]+f[11]+f[12])
            +(T)8*(f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18]);
        moments[2] = 12*f[0]-4*f[1]-4*f[2]-4*f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]-4*f[10]-4*f[11]-4*f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
        moments[3] = -f[1]-f[4]-f[5]-f[6]-f[7]+f[10]+f[13]+f[14]+f[15]+f[16];
        moments[4] = 4*f[1]-f[4]-f[5]-f[6]-f[7]-4*f[10]+f[13]+f[14]+f[15]+f[16];
        moments[5] = -f[2]-f[4]+f[5]-f[8]-f[9]+f[11]+f[13]-f[14]+f[17]+f[18];
        moments[6] = 4*f[2]-f[4]+f[5]-f[8]-f[9]-4*f[11]+f[13]-f[14]+f[17]+f[18];
        moments[7] = -f[3]-f[6]+f[7]-f[8]+f[9]+f[12]+f[15]-f[16]+f[17]-f[18];
        moments[8] = 4*f[3]-f[6]+f[7]-f[8]+f[9]-4*f[12]+f[15]-f[16]+f[17]-f[18];
        moments[9] = 2*f[1]-f[2]-f[3]+f[4]+f[5]+f[6]+f[7]-2*f[8]-2*f[9]+2*f[10]-f[11]-f[12]+f[13]+f[14]+f[15]+f[16]-2*f[17]-2*f[18];
        moments[10] = -4*f[1]+2*f[2]+2*f[3]+f[4]+f[5]+f[6]+f[7]-2*f[8]-2*f[9]-4*f[10]+2*f[11]+2*f[12]+f[13]+f[14]+f[15]+f[16]-2*f[17]-2*f[18];
        moments[11] = f[2]-f[3]+f[4]+f[5]-f[6]-f[7]+f[11]-f[12]+f[13]+f[14]-f[15]-f[16];
        moments[12] = -2*f[2]+2*f[3]+f[4]+f[5]-f[6]-f[7]-2*f[11]+2*f[12]+f[13]+f[14]-f[15]-f[16];
        moments[13] = f[4]-f[5]+f[13]-f[14];
        moments[14] = f[8]-f[9]+f[17]-f[18];
        moments[15] = f[6]-f[7]+f[15]-f[16];
        moments[16] = -f[4]-f[5]+f[6]+f[7]+f[13]+f[14]-f[15]-f[16];
        moments[17] = f[4]-f[5]-f[8]-f[9]-f[13]+f[14]+f[17]+f[18];
        moments[18] = -f[6]+f[7]+f[8]-f[9]+f[15]-f[16]-f[17]+f[18];
        
    }
    
    static void computef_InvM_Smoments(Array<T,19>& f, const Array<T,19> &moments, const T &omega) 
    {
        T mom0 = moments[0] * MRTDescriptor::S[0];
        T mom1 = moments[1] * MRTDescriptor::S[1];
        T mom2 = moments[2] * MRTDescriptor::S[2];
        T mom3 = moments[3] * MRTDescriptor::S[3];
        T mom4 = moments[4] * MRTDescriptor::S[4];
        T mom5 = moments[5] * MRTDescriptor::S[5];
        T mom6 = moments[6] * MRTDescriptor::S[6];
        T mom7 = moments[7] * MRTDescriptor::S[7];
        T mom8 = moments[8] * MRTDescriptor::S[8];
        T mom9 = moments[9] * omega;
        T mom10 = moments[10] * MRTDescriptor::S[10];
        T mom11 = moments[11] * omega;
        T mom12 = moments[12] * MRTDescriptor::S[12];
        T mom13 = moments[13] * omega;
        T mom14 = moments[14] * omega;
        T mom15 = moments[15] * omega;
        T mom16 = moments[16] * MRTDescriptor::S[16];
        T mom17 = moments[17] * MRTDescriptor::S[17];
        T mom18 = moments[18] * MRTDescriptor::S[18];
        
        T mom0tmp = mom0 / (T)19;
        
        f[0] -= mom0tmp-5/(T)399*mom1+1/(T)21*mom2;
        
        T mom1tmp = (T)11/(T)2394*mom1;
        T mom2tmp = mom2/(T)63;
        T mom3_m4 = (T)0.1*(mom3-mom4);
        T mom9_m10 = (mom9-mom10)/(T)18;
        f[1] -= mom0tmp-mom1tmp-mom2tmp-mom3_m4+mom9_m10;
        f[10] -= mom0tmp-mom1tmp-mom2tmp+mom3_m4+mom9_m10;
        
        T mom5_m6 = (T)0.1*(mom5-mom6);
        T mom7_m8 = (T)0.1*(mom7-mom8);
        mom9_m10 *= (T)0.5;
        T mom11_m12 = (mom11-mom12)/(T)12;
        f[2] -= mom0tmp-mom1tmp-mom2tmp-mom5_m6-mom9_m10+mom11_m12;
        f[3] -= mom0tmp-mom1tmp-mom2tmp-mom7_m8-mom9_m10-mom11_m12;
        f[11] -= mom0tmp-mom1tmp-mom2tmp+mom5_m6-mom9_m10+mom11_m12;
        f[12] -= mom0tmp-mom1tmp-mom2tmp+mom7_m8-mom9_m10-mom11_m12;
        
        mom1tmp = (T)4/(T)1197*mom1;
        mom2tmp *= (T)0.25;
        T mom5_p6 = (T)0.1*mom5+(T)0.025*mom6;
        T mom7_p8 = (T)0.1*mom7+(T)0.025*mom8;
        T mom9_p10 = (mom9+mom10*(T)0.5)/(T)18;
        mom14 *= (T)0.25;
        mom17 *= (T)0.125;
        mom18 *= (T)0.125;
        f[8] -= mom0tmp+mom1tmp+mom2tmp-mom5_p6-mom7_p8-mom9_p10+mom14-mom17+mom18;
        f[9] -= mom0tmp+mom1tmp+mom2tmp-mom5_p6+mom7_p8-mom9_p10-mom14-mom17-mom18;
        f[17] -= mom0tmp+mom1tmp+mom2tmp+mom5_p6+mom7_p8-mom9_p10+mom14+mom17-mom18;
        f[18] -= mom0tmp+mom1tmp+mom2tmp+mom5_p6-mom7_p8-mom9_p10-mom14+mom17+mom18;
        
        T mom3_p4 = (T)0.1*mom3+(T)0.025*mom4;
        mom9_p10 *= (T)0.5;
        T mom11_p12 = (mom11+(T)0.5*mom12)/(T)12;
        mom13 *= (T)0.25;
        mom16 *= (T)0.125;
        f[4] -= mom0tmp+mom1tmp+mom2tmp-mom3_p4-mom5_p6+mom9_p10+mom11_p12+mom13-mom16+mom17;
        f[5] -= mom0tmp+mom1tmp+mom2tmp-mom3_p4+mom5_p6+mom9_p10+mom11_p12-mom13-mom16-mom17;
        f[13] -= mom0tmp+mom1tmp+mom2tmp+mom3_p4+mom5_p6+mom9_p10+mom11_p12+mom13+mom16-mom17;
        f[14] -= mom0tmp+mom1tmp+mom2tmp+mom3_p4-mom5_p6+mom9_p10+mom11_p12-mom13+mom16+mom17;
        
        mom15 *= (T)0.25;
        f[15] -= mom0tmp+mom1tmp+mom2tmp+mom3_p4+mom7_p8+mom9_p10-mom11_p12+mom15-mom16+mom18;
        f[16] -= mom0tmp+mom1tmp+mom2tmp+mom3_p4-mom7_p8+mom9_p10-mom11_p12-mom15-mom16-mom18;
        f[6] -= mom0tmp+mom1tmp+mom2tmp-mom3_p4-mom7_p8+mom9_p10-mom11_p12+mom15+mom16-mom18;
        f[7] -= mom0tmp+mom1tmp+mom2tmp-mom3_p4+mom7_p8+mom9_p10-mom11_p12-mom15+mom16+mom18;
    }
    
    static void computeMneqInPlace(Array<T,19> &moments, const Array<T,19> &momentsEq) {
        moments[0] -= momentsEq[0];
        moments[1] -= momentsEq[1];
        moments[2] -= momentsEq[2];
        moments[3] -= momentsEq[3];
        moments[4] -= momentsEq[4];
        moments[5] -= momentsEq[5];
        moments[6] -= momentsEq[6];
        moments[7] -= momentsEq[7];
        moments[8] -= momentsEq[8];
        moments[9] -= momentsEq[9];
        moments[10] -= momentsEq[10];
        moments[11] -= momentsEq[11];
        moments[12] -= momentsEq[12];
        moments[13] -= momentsEq[13];
        moments[14] -= momentsEq[14];
        moments[15] -= momentsEq[15];
        moments[16] -= momentsEq[16];
        moments[17] -= momentsEq[17];
        moments[18] -= momentsEq[18];
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f, const T &omega )
    {
        Array<T,19> moments, momentsEq;
        
        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,3> j(moments[MRTDescriptor::momentumIndexes[0]],moments[MRTDescriptor::momentumIndexes[1]],moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
        
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,3> & j,
                           const T omega )
    {
        Array<T,19> moments, momentsEq;
        
        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        computeEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    /// MRT collision step
    static T smagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,3> & j,
                           const T &omega, const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;
        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeSmagorinskyEquilibriumMoments( momentsEq, rhoBar, j, jSqr, strain, cSmago );
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    /// MRT collision step
    static T smagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                           const T &omega, const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;
        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,3> j(moments[MRTDescriptor::momentumIndexes[0]],moments[MRTDescriptor::momentumIndexes[1]],moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeSmagorinskyEquilibriumMoments( momentsEq, rhoBar, j, jSqr, strain, cSmago );
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    static void addGuoForce( Array<T,Descriptor::q>& f, const Array<T,Descriptor::d>& g,
                             Array<T,Descriptor::d> const& u,
                             const T &omega, T amplitude )
    {
        Array<T,Descriptor::q> forcing, momForce;
        T gx = amplitude * g[0];
        T gy = amplitude * g[1];
        T gz = amplitude * g[2];

        T g_u   = gx * u[0] + gy * u[1] + gz * u[2];
        
        momForce[0] = 0;
        momForce[1] = 38*g_u;
        momForce[2] = -11*g_u;
        momForce[3] = gx;
        momForce[4] = -(T)2/3*gx;
        momForce[5] = gy;
        momForce[6] = -(T)2/3*gy;
        momForce[7] = gz;
        momForce[8] = -(T)2/3*gz;
        momForce[9] = 4*gx*u[0]-2*gy*u[1]-2*gz*u[2];
        momForce[10] = -2*gx*u[0]+gy*u[1]+gz*u[2];
        momForce[11] = 2*gy*u[1]-2*gz*u[2];
        momForce[12] = -gy*u[1]+gz*u[2];
        momForce[13] = gx*u[1]+gy*u[0];
        momForce[14] = gy*u[2]+gz*u[1];
        momForce[15] = gx*u[2]+gz*u[0];
        momForce[16] = 0;
        momForce[17] = 0;
        momForce[18] = 0;
        
        momForce[0] *= (T)0.5;
        momForce[1] *= (T)0.5;
        momForce[2] *= (T)0.5;
        momForce[3] *= (T)0.5;
        momForce[4] *= (T)0.5;
        momForce[5] *= (T)0.5;
        momForce[6] *= (T)0.5;
        momForce[7] *= (T)0.5;
        momForce[8] *= (T)0.5;
        momForce[9] *= (T)0.5;
        momForce[10] *= (T)0.5;
        momForce[11] *= (T)0.5;
        momForce[12] *= (T)0.5;
        momForce[13] *= (T)0.5;
        momForce[14] *= (T)0.5;
        momForce[15] *= (T)0.5;
        momForce[16] *= (T)0.5;
        momForce[17] *= (T)0.5;
        momForce[18] *= (T)0.5;
        
        computef_InvM_Smoments(f, momForce, omega);
        
        static const T oneOver6 = (T)1/(T)6;
        static const T oneOver12 = (T)1/(T)12;
        
        f[0]  += -g_u;
        
        f[1]  +=  oneOver6*(gx*(-(T)1+2*u[0])-gy*u[1]-gz*u[2]);
        f[10] +=  oneOver6*(gx*((T)1+2*u[0]) -gy*u[1] -gz*u[2]);
        
        f[2]  += -oneOver6*(gx*u[0]+gy*((T)1-2*u[1])+gz*u[2]);
        f[3]  += -oneOver6*(gx*u[0] + gy*u[1] + gz*((T)1-2*u[2]));
        f[11] += -oneOver6*(gx*u[0] +gy*(-(T)1-2*u[1]) +gz*u[2]);
        f[12] += -oneOver6*(gx*u[0] + gy*u[1] + gz*(-(T)1-2*u[2]));
        
        f[4]  +=  oneOver12*( gx*(-(T)1+2*u[0]+3*u[1]) + gy*(-(T)1+2*u[1]+3*u[0]) - gz*u[2]);
        f[5]  +=  oneOver12*( gx*(-(T)1+2*u[0]-3*u[1]) + gy*((T)1+2*u[1]-3*u[0]) - gz*u[2]);
        f[6]  +=  oneOver12*(gx*(-(T)1+2*u[0]+3*u[2]) - gy*u[1] + gz*(-(T)1+2*u[2]+3*u[0]));
        f[7]  +=  oneOver12*(gx*(-(T)1+2*u[0]-3*u[2]) - gy*u[1] + gz*((T)1+2*u[2]-3*u[0]));
        f[8]  += -oneOver12*(gx*u[0] + gy*((T)1-2*u[1]-3*u[2]) + gz*((T)1-2*u[2]-3*u[1]));
        f[9]  += -oneOver12*(gx*u[0] + gy*((T)1-2*u[1]+3*u[2]) + gz*(-(T)1-2*u[2]+3*u[1]));
        f[13] +=  oneOver12*(gx*((T)1+2*u[0]+3*u[1]) + gy*((T)1+2*u[1]+3*u[0]) - gz*u[2]);
        f[14] +=  oneOver12*(gx*((T)1+2*u[0]-3*u[1]) + gy*(-(T)1+2*u[1]-3*u[0]) - gz*u[2]);
        f[15] +=  oneOver12*(gx*((T)1+2*u[0]+3*u[2]) - gy*u[1] + gz*((T)1+2*u[2]+3*u[0]));
        f[16] +=  oneOver12*(gx*((T)1+2*u[0]-3*u[2]) - gy*u[1] + gz*(-(T)1+2*u[2]-3*u[0]));
        f[17] += -oneOver12*(gx*u[0] + gy*(-(T)1-2*u[1]-3*u[2]) + gz*(-(T)1-2*u[2]-3*u[1]));
        f[18] += -oneOver12*(gx*u[0] + gy*(-(T)1-2*u[1]+3*u[2]) + gz*((T)1-2*u[2]+3*u[1]));
    }
    
    /// MRT collision step
    static T mrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    const T &omega, 
                                    const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = mrtCollision( f, rhoBar, j, omega );
        addGuoForce( f, force, u, omega, amplitude );
        
        return jSqr;
    }
    
    /// Smagorinsky MRT collision step
    static T smagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                                       const T &rhoBar, const Array<T,Descriptor::d> & u,
                                                       T invM_S[Descriptor::q][Descriptor::q], 
                                                       const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n > &strain, T cSmago, 
                                                       const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = smagorinskyMrtCollision( f, rhoBar, j, invM_S, strain, cSmago );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    /// MRT collision step
    static T incMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    T invM_S[Descriptor::q][Descriptor::q], 
                                    const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = incMrtCollision( f, rhoBar, j, invM_S );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    /// Smagorinsky MRT collision step
    static T incSmagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                                       const T &rhoBar, const Array<T,Descriptor::d> & u,
                                                       T invM_S[Descriptor::q][Descriptor::q], 
                                                       const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n > &strain, T cSmago, 
                                                       const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = incSmagorinskyMrtCollision( f, rhoBar, j, invM_S, strain, cSmago );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncEquilibriumMoments( Array<T,Descriptor::q>& momentsEq,
                                            T rhoBar, Array<T,3> const& j, T jSqr )
    {
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*jSqr-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*jSqr+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*j[0]*j[0]-j[1]*j[1]-j[2]*j[2]);
        momentsEq[10] = (-j[0]*j[0]+(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2]);
        momentsEq[11] = (j[1]*j[1]-j[2]*j[2]);
        momentsEq[12] = (-(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2]);
        momentsEq[13] = j[1]*j[0];
        momentsEq[14] = j[2]*j[1];
        momentsEq[15] = j[2]*j[0];
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeIncSmagorinskyEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                                       T rhoBar, Array<T,3> const& j, T jSqr, const Array<T,6> &strain, T cSmago )
    {
        typedef SymmetricTensorImpl<T,3> S;
        T sNorm = sqrt((T)2*SymmetricTensorImpl<T,3>::tensorNormSqr(strain));
        T smagoFactor = (T)2*cSmago * cSmago * sNorm;
        
        T ux2 = smagoFactor * strain[S::xx];
        T uy2 = smagoFactor * strain[S::yy];
        T uz2 = smagoFactor * strain[S::zz];
        
        T uxuy = smagoFactor * strain[S::xy];
        T uyuz = smagoFactor * strain[S::yz];
        T uxuz = smagoFactor * strain[S::xz];
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*(jSqr + ux2+uy2+uz2) -(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*(jSqr + ux2+uy2+uz2)+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*(j[0]*j[0]+ux2)-(j[1]*j[1]+uy2)-(j[2]*j[2]+uz2));
        momentsEq[10] = (-(j[0]*j[0]+ux2)+(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2));
        momentsEq[11] = (j[1]*j[1]+uy2-(j[2]*j[2]+uz2));
        momentsEq[12] = (-(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2));
        momentsEq[13] = j[1]*j[0]+uxuy;
        momentsEq[14] = j[2]*j[1]+uyuz;
        momentsEq[15] = j[2]*j[0]+uxuz;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// MRT collision step
    static T incMrtCollision( Array<T,Descriptor::q>& f,
                                   const T &omega )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,3> j(moments[MRTDescriptor::momentumIndexes[0]],moments[MRTDescriptor::momentumIndexes[1]],moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        computeIncEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    /// MRT collision step
    static T incMrtCollision( Array<T,Descriptor::q>& f,
                                   const T &rhoBar, const Array<T,3> & j,
                                   const T &omega )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        computeIncEquilibriumMoments(momentsEq,rhoBar,j,jSqr);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    /// MRT collision step
    static T incSmagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                                              const T &omega, const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;
        
        computeMoments(moments,f);
        T rhoBar = moments[0];
        Array<T,3> j(moments[MRTDescriptor::momentumIndexes[0]],moments[MRTDescriptor::momentumIndexes[1]],moments[MRTDescriptor::momentumIndexes[2]]);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeIncSmagorinskyEquilibrium(momentsEq,rhoBar,j,jSqr,strain,cSmago);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }
    
    /// MRT collision step
    static T incSmagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                                              const T &rhoBar, const Array<T,3> & j,
                                              const T &omega, const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeIncSmagorinskyEquilibrium(momentsEq,rhoBar,j,jSqr,strain,cSmago);
        computeMneqInPlace(moments,momentsEq); // moments become mNeq
        computef_InvM_Smoments(f, moments, omega);
        
        return jSqr;
    }

};


}  // namespace plb

#endif
