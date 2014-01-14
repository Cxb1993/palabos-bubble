/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2014 FlowKit Sarl
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


#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cstdio>
#include <cmath>

#undef POISEUILLE

using namespace plb;

typedef double T;

#define PADDING 8

//#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
#define DESCRIPTOR descriptors::ForcedD3Q27Descriptor


T fltEps = 1.0e-5;
T dblEps = 1.0e-9;

std::string outDir("tmp/");


// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param {
    T radius;                       // Radius of the orifice from which the water jet emerges.
    T radius_LB;
    T lx, ly, lz;                   // Dimensions of the simulation domain.
    T orifice_height;               // Height of the orifice.
    T orifice_height_LB;
    T orifice_thickness;
    T orifice_thickness_LB;
    T initial_fluid_height;         // Initial "thickness" of the water inside the orifice.
    T initial_fluid_height_LB;
    T water_pool_height;            // Initial "thickness" of the water at the bottom of the simulation domain.
    T water_pool_height_LB;

    bool useBubblePressureModel;    // Use a model where the bubble pressure depends on the bubble volume.
    T bubbleVolumeRatio;            // Parameter of the bubble-pressure model. Must be >= 1.0.
    T alpha, beta;                  // Parameters of the bubble-pressure model.

    plint resolution;               // Number of nodes in the radius of the orifice.
    T u_LB;                         // u_LB = U_inlet * dt / dx. From this, the time step is calculated.
    plint nx, ny, nz;
    plint maxIter;                  // Total number of iterations.
    plint statIter;                 // Number of iterations for statistics output.
    plint outIter;                  // Number of iterations for disk output.

    T rho;                          // Density of water.
    T rho_LB;

    T Re;                           // Reynolds number: Re = U_inlet * R / nu.
    T We;                           // Weber number: We = rho * U_inlet^2 * R / gamma. "gamma" is the surface tension.
    T Bo;                           // Bond number: Bo = rho * g * R^2 / gamma.

    T contact_angle;                // In degrees, if < 0 then the contact angle algorithm is deactivated.
    Array<T,3> external_force_LB;

    T inlet_velocity;               // Velocity of water at the inlet of the orifice.
    T inlet_velocity_LB;

    T cSmago;                       // Parameter of the Smagorinsky LES model.

    T surface_tension_LB;

    T omega;

    T dx, dt;
} param;


void setPhysicalParam(std::string xmlInputFileName)
{
    XMLreader document(xmlInputFileName);

    document["geometry"]["simulationDomain"]["lx"].read(param.lx);
    document["geometry"]["simulationDomain"]["ly"].read(param.ly);
    document["geometry"]["simulationDomain"]["lz"].read(param.lz);
    document["geometry"]["radius"].read(param.radius);
    document["geometry"]["orificeHeight"].read(param.orifice_height);
    document["geometry"]["orificeThickness"].read(param.orifice_thickness);
    document["geometry"]["initialFluidHeight"].read(param.initial_fluid_height);
    document["geometry"]["waterPoolHeight"].read(param.water_pool_height);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["Re"].read(param.Re);
    document["fluid"]["We"].read(param.We);
    document["fluid"]["Bo"].read(param.Bo);
    document["fluid"]["contactAngle"].read(param.contact_angle);

    document["numerics"]["resolution"].read(param.resolution);
    document["numerics"]["uLB"].read(param.u_LB);
    document["numerics"]["inletVelocity"].read(param.inlet_velocity);
    document["numerics"]["maxIter"].read(param.maxIter);
    document["numerics"]["cSmago"].read(param.cSmago);
    document["numerics"]["useBubblePressureModel"].read(param.useBubblePressureModel);
    document["numerics"]["bubbleVolumeRatio"].read(param.bubbleVolumeRatio);
    document["numerics"]["alpha"].read(param.alpha);
    document["numerics"]["beta"].read(param.beta);

    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);

    pcout << "radius = " << param.radius << std::endl;
    pcout << "lx = " << param.lx << std::endl;
    pcout << "ly = " << param.ly << std::endl;
    pcout << "lz = " << param.lz << std::endl;
    pcout << "orifice_height = " << param.orifice_height << std::endl;
    pcout << "orifice_thickness = " << param.orifice_thickness << std::endl;
    pcout << "initial_fluid_height = " << param.initial_fluid_height << std::endl;
    pcout << "water_pool_height = " << param.water_pool_height << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;
    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "u_LB = " << param.u_LB << std::endl;
    pcout << "rho = " << param.rho << std::endl;
    pcout << "Re = " << param.Re << std::endl;
    pcout << "We = " << param.We << std::endl;
    pcout << "Bo = " << param.Bo << std::endl;
    pcout << "contact_angle = " << param.contact_angle << std::endl;
    pcout << "inlet_velocity = " << param.inlet_velocity << std::endl;
    pcout << "cSmago = " << param.cSmago << std::endl;
    pcout << "useBubblePressureModel = " << (param.useBubblePressureModel ? "true" : "false") << std::endl;
    if (param.useBubblePressureModel) {
        pcout << "bubbleVolumeRatio = " << param.bubbleVolumeRatio << std::endl;
        pcout << "alpha = " << param.alpha << std::endl;
        pcout << "beta = " << param.beta << std::endl;
    }
    pcout << std::endl;
}


void computeLBParam(void)
{
    // Derived quantities

    param.dx = param.radius / (param.resolution - 1.0);
    param.dt = (param.u_LB / param.inlet_velocity) * param.dx;

    param.radius_LB = util::roundToInt(param.radius / param.dx);

    param.nx = util::roundToInt(param.lx / param.dx);
    param.ny = util::roundToInt(param.ly / param.dx);
    param.nz = util::roundToInt(param.lz / param.dx);

    param.orifice_height_LB = util::roundToInt(param.orifice_height / param.dx);
    param.orifice_thickness_LB = util::roundToInt(param.orifice_thickness / param.dx);
    param.initial_fluid_height_LB = util::roundToInt(param.initial_fluid_height / param.dx);
    param.water_pool_height_LB = util::roundToInt(param.water_pool_height / param.dx);
    PLB_ASSERT(param.orifice_height_LB > 0 && fabs(param.orifice_height_LB) > fltEps)
    PLB_ASSERT(param.orifice_thickness_LB > 0 && fabs(param.orifice_thickness_LB) > fltEps)
    PLB_ASSERT(param.orifice_height_LB - param.initial_fluid_height_LB >= 0 ||
            fabs(param.orifice_height_LB - param.initial_fluid_height_LB) <= fltEps)
    PLB_ASSERT(param.water_pool_height_LB > 0 && fabs(param.water_pool_height_LB) > fltEps)

    param.inlet_velocity_LB = param.inlet_velocity * param.dt / param.dx;

    param.rho_LB = 1.0;

    T nu_LB = param.inlet_velocity_LB * param.radius_LB / param.Re;
    param.omega = 1.0 / (3.0 * nu_LB + 0.5);

    param.surface_tension_LB = param.rho_LB * param.inlet_velocity_LB * param.inlet_velocity_LB * param.radius_LB / param.We;

    T g_LB = param.surface_tension_LB * param.Bo / (param.rho_LB * param.radius_LB * param.radius_LB);
    param.external_force_LB = Array<T,3>((T) 0, (T) 0, -g_LB);
}


int initialFlags(plint iX, plint iY, plint iZ)
{
    if (iX == 0 || iX == param.nx - 1 || iY == 0 || iY == param.ny - 1 || iZ == 0 || iZ == param.nz - 1) {
        return twoPhaseFlag::wall;
    }

    T cx = 0.5 * param.nx;
    T cy = 0.5 * param.ny;

    T dx = iX - cx;
    T dy = iY - cy;

    T rSmallSqr = param.radius_LB * param.radius_LB;
    T rBigSqr = (param.radius_LB + param.orifice_thickness_LB) * (param.radius_LB + param.orifice_thickness_LB);
    T rSqr = dx * dx + dy * dy;

    int inside_wall = (iZ >= param.nz - param.orifice_height_LB && rSqr > rSmallSqr && rSqr <= rBigSqr);
    inside_wall = inside_wall || (iZ == param.nz - 2 && rSqr > rSmallSqr); // Second layer of wall nodes for inlet BC.
    int inside_fluid = iZ >= param.nz - param.initial_fluid_height_LB && rSqr <= rSmallSqr;
    inside_fluid = inside_fluid || iZ <= param.water_pool_height_LB;

    if (inside_wall) {
        return twoPhaseFlag::wall;
    } else if (inside_fluid) {
        return twoPhaseFlag::fluid;
    } else {
        return twoPhaseFlag::empty;
    }
}


#ifdef POISEUILLE
class InletVelocity {
public:
    void operator()(T x, T y, T z, Array<T,3>& u) const
    {
        T xc = 0.5 * param.nx;
        T yc = 0.5 * param.ny;

        T xp = x - xc;
        T yp = y - yc;

        T r = sqrt(xp * xp + yp * yp);

        if (r <= param.radius_LB) {
            u[0] =  0.0;
            u[1] =  0.0;
            u[2] = -2.0 * param.inlet_velocity_LB * (1.0 - (r / param.radius_LB) * (r / param.radius_LB));
        } else {
            u.resetToZero();
        }
    }
};
#endif


void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR> *fields, MultiScalarField3D<plint> *tagMatrix, plint iter)
{
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox()));

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox().enlarge(0));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "smoothedInterface_", iter, PADDING)+".stl");
    }
    triangles.clear();
    isoSurfaceMarchingCube(triangles, fields->volumeFraction, isoLevels, fields->volumeFraction.getBoundingBox().enlarge(0));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "interface_", iter, PADDING)+".stl");
    }

    //T coef = 1.0 / 3.0;
    VtkImageOutput3D<T> vtkOut(createFileName(outDir + "volumeData_", iter, PADDING), param.dx);
    std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(fields->lattice);
    std::auto_ptr<MultiScalarField3D<T> > rho = computeDensity(fields->lattice);
    vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
    //vtkOut.writeData<float>(*rho, "pressure", param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
    vtkOut.writeData<float>(fields->volumeFraction, "volumeFraction", 1.0);
    vtkOut.writeData<float>(*smoothVF, "smoothedVolumeFraction", 1.0);
    vtkOut.writeData<float>(*copyConvert<int,double>(fields->flag, fields->flag.getBoundingBox()), "flag", 1.0);
    //vtkOut.writeData<float>(fields->outsideDensity, "outsidePressure",
    //        param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
    if (tagMatrix != 0) {
        vtkOut.writeData<float>(*copyConvert<plint,double>(*tagMatrix), "bubbleTags", 1.0);
    }
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    std::cout.precision(10);
    std::scientific(std::cout);

    // Command-line arguments

    if (argc != 2) {
        pcout << "Usage: " << argv[0] << " xml-input-file-name" << std::endl;
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);

    // Set the simulation parameters.

    setPhysicalParam(xmlInputFileName);

    computeLBParam();

    pcout << "Tau: " << 1.0 / param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "nx = " << param.nx << std::endl;
    pcout << "ny = " << param.ny << std::endl;
    pcout << "nz = " << param.nz << std::endl;
    pcout << "Radius (lu) = " << param.radius_LB << std::endl;
    pcout << "Inlet velocity (lu) = " << param.inlet_velocity_LB << std::endl;
    pcout << "Gravity (lu) = " << param.external_force_LB[2] << std::endl;
    pcout << "Surface tension (lu) = " << param.surface_tension_LB << std::endl;

    // Free surface fields construction.

    SparseBlockStructure3D blockStructure(createRegularDistribution3D(param.nx, param.ny, param.nz));

    Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
    pcout << "Dynamics model: Smagorinsky BGK" << std::endl;

    bool useExtendedEnvelopes = false;

    FreeSurfaceFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), param.rho_LB,
            param.surface_tension_LB, param.contact_angle, param.external_force_LB,
            useExtendedEnvelopes);

    fields.periodicityToggleAll(false);

    setToFunction(fields.flag, fields.flag.getBoundingBox(), initialFlags);

    fields.defaultInitialize();

    Box3D inlet;
    inlet.x0 = 1;
    inlet.x1 = param.nx - 2;
    inlet.y0 = 1;
    inlet.y1 = param.ny - 2;
    inlet.z0 = param.nz - 2;
    inlet.z1 = param.nz - 2;

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    boundaryCondition->addVelocityBoundary2P(inlet, fields.lattice);

#ifdef POISEUILLE
    setBoundaryVelocity(fields.lattice, inlet, InletVelocity());
#else
    setBoundaryVelocity(fields.lattice, inlet, Array<T,3>((T) 0, (T) 0, -param.inlet_velocity_LB));
#endif

    delete boundaryCondition;

    defineDynamics(fields.lattice, fields.flag, inlet, new BounceBack<T,DESCRIPTOR>(param.rho_LB),
            (int) twoPhaseFlag::wall);

    BubbleMatch3D *bubbleMatch = 0;
    BubbleHistory3D<T> *bubbleHistory = 0;
    FILE *fp = 0;

    if (param.useBubblePressureModel) {
        bubbleMatch = new BubbleMatch3D(fields.flag);
        bubbleHistory = new BubbleHistory3D<T>(fields.flag);
        std::string fname = outDir + "bubbles.log";
        fp = fopen(fname.c_str(), "w");
        PLB_ASSERT(fp != 0);
    }

    pcout << std::endl;

    pcout << "Starting simulation." << std::endl;
    for (plint iter = 0; iter < param.maxIter; iter++) {
        if ((iter % param.statIter == 0 && iter != 0) || iter == param.maxIter - 1) {
            pcout << "At iteration " << iter << ", t = " << iter * param.dt << std::endl;
            T avE = computeAverageEnergy(fields.lattice);
            avE *= param.dx * param.dx / (param.dt * param.dt);
            pcout << "Average kinetic energy: " << avE << std::endl;
            plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
            pcout << "Number of interface cells: " << numIntCells << std::endl;
            pcout << std::endl;
        }

        if (iter % param.outIter == 0 || iter == param.maxIter - 1) {
            pcout << "Output to disk at iteration " << iter << " ... ";
            if (param.useBubblePressureModel) {
                writeResults(&fields, bubbleMatch->getTagMatrix(), iter);
            } else {
                writeResults(&fields, 0, iter);
            }
            pcout << "done." << std::endl;
            pcout << std::endl;
        }

        // Execute all the free surface data processors.
        if (param.useBubblePressureModel) {
            T bubbleVolumeRatio = iter == 0 ? 1.0 : param.bubbleVolumeRatio;
            bubbleMatch->execute(fields.flag, fields.volumeFraction);
            bubbleHistory->transition(*bubbleMatch, iter, bubbleVolumeRatio);
            bubbleHistory->updateBubblePressure(fields.outsideDensity, param.rho_LB, param.alpha, param.beta);
        }

        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        applyProcessingFunctional(new InletConstVolumeFraction3D<T,DESCRIPTOR>(1.002), inlet, fields.twoPhaseArgs);

        if (param.useBubblePressureModel) {
            if (iter == 0) {
                bubbleHistory->freeze();
            }
            if (iter % param.statIter == 0 || iter == param.maxIter-1) {
                bubbleHistory->timeHistoryLog(outDir + "/bubbleTimeHistory.log");
                bubbleHistory->fullBubbleLog(outDir + "/fullBubbleRecord.log");

                // We do not log frozen bubbles.
                std::map<plint,BubbleInfo3D>::const_iterator it = bubbleHistory->getBubbles().begin();
                T totalBubbleVolume = T();
                plint numBubbles = 0;
                TriangleSet<T> bubbles;
                T pi = acos((T) -1);
                for (; it != bubbleHistory->getBubbles().end(); ++it) {
                    if (it->second.isFrozen()) {
                        continue;
                    }
                    numBubbles++;
                    T v = it->second.getVolume() * param.dx * param.dx * param.dx;
                    T r = pow(v * 3.0 / 4.0 / pi, 1.0 / 3.0);
                    Array<T,3> center = it->second.getCenter() * param.dx;
                    pcout << "Bubble " << it->first
                          << " has volume " << v
                          << ", radius " << r
                          << " and center (" << center[0] << "," << center[1] << "," << center[2] << ")" << std::endl;
                    totalBubbleVolume += v;
                    fprintf(fp, "Time: % .8e, Bubble id: %8ld, Volume: % .8e, Radius: % .8e, Center: (% .8e, % .8e, % .8e)\n",
                            iter * param.dt, it->first, v, r, center[0], center[1], center[2]);

                    // Save STL files of bubbles.

                    if (!util::fpequal(r, (T) 0.0, std::numeric_limits<T>::epsilon())) {
                        plint minNumOfTriangles = 64;
                        TriangleSet<T> triangleSet = constructSphere(center, r, minNumOfTriangles);
                        bubbles.append(triangleSet);
                    }
                }
                std::string fname = createFileName(outDir + "/bubbles_", iter, PADDING) + ".stl";
                bubbles.writeBinarySTL(fname);
                pcout << "At iteration " << iter << ", the number of bubbles is " << numBubbles << std::endl;
                pcout << "The total volume of bubbles is: " << totalBubbleVolume << std::endl;
                pcout << std::endl;
                fflush(fp);
            }
        }
    }

    if (param.useBubblePressureModel) {
        delete bubbleHistory;
        delete bubbleMatch;
        fclose(fp);
    }

    exit(0);
}

