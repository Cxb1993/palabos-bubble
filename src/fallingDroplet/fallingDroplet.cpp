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
#include <vector>
#include <string>

using namespace plb;

#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
//#define DESCRIPTOR descriptors::ForcedD3Q27Descriptor

#define PADDING 8

typedef double T;

std::string outDir("./tmp/");

struct SimulationParameters {
    /*
     * Parameters set by the user.
     */

    T lx, ly, lz;                                   // Geometric parameters.
    T fluidPoolHeight;
    Array<T,3> sphereCenter;
    Array<T,3> sphereCenter2;
    T sphereRadius;

    T rho;                                          // Material parameters.
    T Re;
    T Ga;
    T La;
    T contactAngle;

    plint resolution;                               // Numerical parameters.
    T u_LB;
    plint maxIter;
    T cSmago;
    bool useBubblePressureModel;
    bool freezeLargestBubble;
    T bubbleVolumeRatio;
    T alpha, beta;

    plint statIter;                                 // Output parameters.
    plint outIter;
    bool writeVTK;

    /*
     * Parameters NOT set by the user.
     */

    T dx;
    T dt;
    plint nx, ny, nz;
    Array<T,3> sphereCenter_LB;
    Array<T,3> sphereCenter_LB2;
    T sphereRadius_LB;
    T fluidPoolHeight_LB;
    Array<T,3> initialVelocity_LB;
    Array<T,3> g_LB;
    T rho_LB;
    T nu_LB;
    T sigma_LB;
    T omega;
};

SimulationParameters param;

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    XMLreader document(xmlInputFileName);

    document["geometry"]["simulationDomain"]["lx"].read(param.lx);
    document["geometry"]["simulationDomain"]["ly"].read(param.ly);
    document["geometry"]["simulationDomain"]["lz"].read(param.lz);

    document["geometry"]["fluidPoolHeight"].read(param.fluidPoolHeight);

    std::vector<T> sphereCenter;
    document["geometry"]["sphere"]["center"].read(sphereCenter);
    PLB_ASSERT(sphereCenter.size() == 3);
    param.sphereCenter[0] = sphereCenter[0];
    param.sphereCenter[1] = sphereCenter[1];
    param.sphereCenter[2] = sphereCenter[2];
    document["geometry"]["sphere"]["radius"].read(param.sphereRadius);

    sphereCenter.clear();
    document["geometry"]["sphere"]["secondCenter"].read(sphereCenter);
    PLB_ASSERT(sphereCenter.size() == 3);
    param.sphereCenter2[0] = sphereCenter[0];
    param.sphereCenter2[1] = sphereCenter[1];
    param.sphereCenter2[2] = sphereCenter[2];

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["Re"].read(param.Re);
    document["fluid"]["Ga"].read(param.Ga);
    document["fluid"]["La"].read(param.La);
    document["fluid"]["contactAngle"].read(param.contactAngle);
    T pi = acos((T) -1);
    param.contactAngle *= pi/180.0;

    document["numerics"]["resolution"].read(param.resolution);
    document["numerics"]["uLB"].read(param.u_LB);
    document["numerics"]["maxIter"].read(param.maxIter);
    document["numerics"]["cSmago"].read(param.cSmago);
    document["numerics"]["useBubblePressureModel"].read(param.useBubblePressureModel);
    document["numerics"]["freezeLargestBubble"].read(param.freezeLargestBubble);
    document["numerics"]["bubbleVolumeRatio"].read(param.bubbleVolumeRatio);
    document["numerics"]["alpha"].read(param.alpha);
    document["numerics"]["beta"].read(param.beta);

    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
    document["output"]["writeVTK"].read(param.writeVTK);
}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{
    // Derived quantities.

    param.dx = 1.0 / (param.resolution - 1.0);
    param.dt = param.u_LB * param.dx;

    T l_LB = param.resolution - 1.0;

    param.nx = util::roundToInt(param.lx / param.dx);
    param.ny = util::roundToInt(param.ly / param.dx);
    param.nz = util::roundToInt(param.lz / param.dx);

    param.sphereCenter_LB = (1.0 / param.dx) * param.sphereCenter;
    param.sphereCenter_LB2 = (1.0 / param.dx) * param.sphereCenter2;
    param.sphereRadius_LB = (1.0 / param.dx) * param.sphereRadius;
    param.fluidPoolHeight_LB = (1.0 / param.dx) * param.fluidPoolHeight;

    param.rho_LB = 1.0;

    if (param.Re > 0.0) {
        param.initialVelocity_LB = Array<T,3>((T) 0, (T) 0, -param.u_LB); // U = 1 in dimensionless units.
        param.nu_LB = param.u_LB * l_LB / param.Re;
        if (param.Ga > 0.0) {
            T norm_g_LB = param.Ga * param.nu_LB * param.nu_LB / (l_LB * l_LB * l_LB);
            param.g_LB = Array<T,3>((T) 0, (T) 0, -norm_g_LB);
        } else {
            param.g_LB.resetToZero();
        }
    } else {
        param.initialVelocity_LB.resetToZero();
        if (param.Ga > 0.0) {
            T norm_g_LB = param.dt * param.dt / param.dx; // g = 1 in dimensionless units.
            param.g_LB = Array<T,3>((T) 0, (T) 0, -norm_g_LB);
            param.nu_LB = sqrt(norm_g_LB * l_LB * l_LB * l_LB / param.Ga);
        } else {
            param.g_LB.resetToZero();
            param.nu_LB = param.dt / (param.dx * param.dx); // nu = 1 in dimensionless units.
        }
    }

    T mu_LB = param.nu_LB * param.rho_LB;
    param.sigma_LB = param.La * mu_LB * mu_LB / (param.rho_LB * l_LB);

    param.omega = 1.0 / (3.0 * param.nu_LB + 0.5);
}

void printSimulationParameters(SimulationParameters const& param)
{
    pcout << "lx = " << param.lx << std::endl;
    pcout << "ly = " << param.ly << std::endl;
    pcout << "lz = " << param.lz << std::endl;
    pcout << "fluidPoolHeight = " << param.fluidPoolHeight << std::endl;
    pcout << "sphereCenter = (" << param.sphereCenter[0] << ", "
                                << param.sphereCenter[1] << ", "
                                << param.sphereCenter[2] << ")" << std::endl;
    pcout << "sphereRadius = " << param.sphereRadius << std::endl;

    pcout << "rho = " << param.rho << std::endl;
    pcout << "Re = " << param.Re << std::endl;
    pcout << "Ga = " << param.Ga << std::endl;
    pcout << "La = " << param.La << std::endl;
    pcout << "contractAngle = " << param.contactAngle * 180.0 / acos((T) -1) << std::endl;

    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "uLB = " << param.u_LB << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "cSmago = " << param.cSmago << std::endl;
    pcout << "useBubblePressureModel = " << (param.useBubblePressureModel ? "true" : "false") << std::endl;
    pcout << "freezeLargestBubble = " << (param.freezeLargestBubble ? "true" : "false") << std::endl;
    pcout << "bubbleVolumeRatio = " << param.bubbleVolumeRatio << std::endl;
    pcout << "alpha = " << param.alpha << std::endl;
    pcout << "beta = " << param.beta << std::endl;

    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;

    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "dt / (dx * dx) = " << param.dt / (param.dx * param.dx) << std::endl;
    pcout << "nx = " << param.nx << std::endl;
    pcout << "ny = " << param.ny << std::endl;
    pcout << "nz = " << param.nz << std::endl;
    pcout << "fluidPoolHeight_LB = " << param.fluidPoolHeight_LB << std::endl;
    pcout << "sphereCenter_LB = (" << param.sphereCenter_LB[0] << ", "
                                   << param.sphereCenter_LB[1] << ", "
                                   << param.sphereCenter_LB[2] << ")" << std::endl;
    pcout << "sphereCenter_LB2 = (" << param.sphereCenter_LB2[0] << ", "
                                   << param.sphereCenter_LB2[1] << ", "
                                   << param.sphereCenter_LB2[2] << ")" << std::endl;
    pcout << "sphereRadius_LB = " << param.sphereRadius_LB << std::endl;
    pcout << "initialVelocity_LB = (" << param.initialVelocity_LB[0] << ", "
                                      << param.initialVelocity_LB[1] << ", "
                                      << param.initialVelocity_LB[2] << ")" << std::endl;
    pcout << "g_LB = (" << param.g_LB[0] << ", " << param.g_LB[1] << ", " << param.g_LB[2] << ")" << std::endl;
    pcout << "rho_LB = " << param.rho_LB << std::endl;
    pcout << "nu_LB = " << param.nu_LB << std::endl;
    pcout << "sigma_LB = " << param.sigma_LB << std::endl;
    pcout << "omega = " << param.omega << std::endl;
    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << std::endl;
}

bool insideSphere(T x, T y, T z)
{
    Array<T,3> pos(x, y, z);

    T r = norm<T,3>(pos - param.sphereCenter_LB);
    T r2 = norm<T,3>(pos - param.sphereCenter_LB2);
    if (r <= param.sphereRadius_LB || r2 <= param.sphereRadius_LB) {
        return true;
    }
    return false;
}

bool insideFluidPool(T x, T y, T z)
{
    if (z <= param.fluidPoolHeight_LB) {
        return true;
    }
    return false;
}

bool insideFluid(T x, T y, T z)
{
    if (insideFluidPool(x, y, z) || insideSphere(x, y, z)) {
        return true;
    }
    return false;
}

// Specifies the initial condition for the fluid (each cell is assigned the
// flag "fluid", "empty", or "wall").
int initialFluidFlags(plint iX, plint iY, plint iZ)
{
    if (insideFluid(iX, iY, iZ)) {
        return twoPhaseFlag::fluid;
    }
    return twoPhaseFlag::empty;
}

void writeResults(FreeSurfaceFields3D<T,DESCRIPTOR> *fields, MultiScalarField3D<plint> *tagMatrix, plint iT, bool writeVTK)
{
    pcout << "Writing results" << std::endl;
    pcout << "\tSmoothing volume fraction..." << std::flush;
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox()));

    pcout << "done" << std::endl;

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    pcout << "\tMarching cubes smooth..." << std::flush;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox().enlarge(-2));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "smoothedInterface_", iT, PADDING)+".stl");
    }
    pcout << "done" << std::endl << "\tMarching cubes fraction..." << std::flush;
    triangles.clear();
    isoSurfaceMarchingCube(triangles, fields->volumeFraction, isoLevels, fields->volumeFraction.getBoundingBox().enlarge(-2));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(param.dx);
        triangleSet.writeBinarySTL(createFileName(outDir + "interface_", iT, PADDING)+".stl");
    }
    pcout << "done" << std::endl;

    if (writeVTK)
    {
        // T coef = 1.0 / 3.0;
        VtkImageOutput3D<T> vtkOut(createFileName(outDir + "volumeData_", iT, PADDING), param.dx);
        pcout << "Computing velocity..." << std::flush;
        std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(fields->lattice);
        pcout << "done" << std::endl << "Computing density..." << std::flush;
        std::auto_ptr<MultiScalarField3D<T> > rho = computeDensity(fields->lattice);
        pcout << "done" << std::endl << "Writing velocity..." << std::flush;
        vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
        //vtkOut.writeData<float>(*rho, "pressure", param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
        pcout << "done" << std::endl << "Writing volume fraction..." << std::flush;
        vtkOut.writeData<float>(fields->volumeFraction, "volumeFraction", 1.0);
        pcout << "done" << std::endl << "Writing smoothed volume fraction..." << std::flush;
        vtkOut.writeData<float>(*smoothVF, "smoothedVolumeFraction", 1.0);
        pcout << "done" << std::endl;
        //vtkOut.writeData<float>(fields->outsideDensity, "outsidePressure",
        //        param.rho * coef * (param.dx * param.dx) / (param.dt * param.dt));
        if (tagMatrix != 0) {
            pcout << "Writing tags..." << std::flush;
            vtkOut.writeData<float>(*copyConvert<plint,double>(*tagMatrix), "bubbleTags", 1.0);
            pcout << "done" << std::endl;
        }
    }
}


int main(int argc, char **argv)
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

    readUserDefinedSimulationParameters(xmlInputFileName, param);
    calculateDerivedSimulationParameters(param);
    printSimulationParameters(param);

    SparseBlockStructure3D blockStructure(createRegularDistribution3D(param.nx, param.ny, param.nz));

    Dynamics<T,DESCRIPTOR>* dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);

    FreeSurfaceFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), param.rho_LB,
            param.sigma_LB, param.contactAngle, param.g_LB, false);

    // Initialization

    pcout << "Setting up initial condition." << std::endl;

    analyticalIniVolumeFraction(fields.volumeFraction, fields.flag, insideFluid, 32);

    Box3D bottom(0, param.nx-1, 0, param.ny-1, 0, 0);
    Box3D top(0, param.nx-1, 0, param.ny-1, param.nz-1, param.nz-1);
    Box3D lateral1(0, 0, 0, param.ny-1, 0, param.nz-1);
    Box3D lateral2(param.nx-1, param.nx-1, 0, param.ny-1, 0, param.nz-1);
    Box3D lateral3(0, param.nx-1, 0, 0, 0, param.nz-1);
    Box3D lateral4(0, param.nx-1, param.ny-1, param.ny-1, 0, param.nz-1);

    setToConstant(fields.flag, bottom, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, top, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral1, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral2, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral3, (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, lateral4, (int) twoPhaseFlag::wall);

    fields.periodicityToggleAll(false);
    fields.partiallyDefaultInitialize();

    plint bx0 = util::roundToInt(param.sphereCenter_LB[0]) - util::roundToInt(param.sphereRadius_LB) - 2;
    plint bx1 = util::roundToInt(param.sphereCenter_LB[0]) + util::roundToInt(param.sphereRadius_LB) + 2;
    plint by0 = util::roundToInt(param.sphereCenter_LB[1]) - util::roundToInt(param.sphereRadius_LB) - 2;
    plint by1 = util::roundToInt(param.sphereCenter_LB[1]) + util::roundToInt(param.sphereRadius_LB) + 2;
    plint bz0 = util::roundToInt(param.sphereCenter_LB[2]) - util::roundToInt(param.sphereRadius_LB) - 2;
    plint bz1 = util::roundToInt(param.sphereCenter_LB[2]) + util::roundToInt(param.sphereRadius_LB) + 2;
    Box3D iniVelBox(bx0, bx1, by0, by1, bz0, bz1);
    applyProcessingFunctional(new ConstantIniVelocityFreeSurface3D<T,DESCRIPTOR>(param.initialVelocity_LB, param.rho_LB),
        iniVelBox, fields.twoPhaseArgs);

    bx0 = util::roundToInt(param.sphereCenter_LB2[0]) - util::roundToInt(param.sphereRadius_LB) - 2;
    bx1 = util::roundToInt(param.sphereCenter_LB2[0]) + util::roundToInt(param.sphereRadius_LB) + 2;
    by0 = util::roundToInt(param.sphereCenter_LB2[1]) - util::roundToInt(param.sphereRadius_LB) - 2;
    by1 = util::roundToInt(param.sphereCenter_LB2[1]) + util::roundToInt(param.sphereRadius_LB) + 2;
    bz0 = util::roundToInt(param.sphereCenter_LB2[2]) - util::roundToInt(param.sphereRadius_LB) - 2;
    bz1 = util::roundToInt(param.sphereCenter_LB2[2]) + util::roundToInt(param.sphereRadius_LB) + 2;
    Box3D iniVelBox2(bx0, bx1, by0, by1, bz0, bz1);
    applyProcessingFunctional(new ConstantIniVelocityFreeSurface3D<T,DESCRIPTOR>(param.initialVelocity_LB, param.rho_LB),
        iniVelBox2, fields.twoPhaseArgs);

    plint iniIter = 0;

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

    // Main iteration loop.

    for (plint iT = iniIter; iT < param.maxIter; iT++) {
        if (iT % param.statIter == 0 || iT == param.maxIter-1) {
            pcout << "At iteration " << iT << ", t = " << iT * param.dt << std::endl;
            T avE = computeAverageEnergy(fields.lattice);
            avE *= param.rho * param.dx * param.dx / (param.dt * param.dt);
            pcout << "Average kinetic energy: " << avE << std::endl;
            plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
            pcout << "Number of interface cells: " << numIntCells << std::endl;
            if (iT != iniIter) {
                pcout << "Time spent for each iteration: "
                      << global::timer("iteration").getTime() / (T) param.statIter << std::endl;
                global::timer("iteration").reset();
            }
            pcout << std::endl;
        }

        if (iT % param.outIter == 0 || iT == param.maxIter-1) {
            pcout << "Writing results at iteration " << iT << ", t = " << iT * param.dt << std::endl;
            global::timer("images").restart();
            if (param.useBubblePressureModel) {
                writeResults(&fields, bubbleMatch->getTagMatrix(), iT, param.writeVTK);
            } else {
                writeResults(&fields, 0, iT, param.writeVTK);
            }
            global::timer("images").stop();
            pcout << "Time spent for writing results: " << global::timer("images").getTime() << std::endl;
            pcout << std::endl;
        }

        global::timer("iteration").start();
        if (param.useBubblePressureModel) {
            T bubbleVolumeRatio = iT == 0 ? 1.0 : param.bubbleVolumeRatio;
            //pcout << "Starting bubble match execute...";
            bubbleMatch->execute(fields.flag, fields.volumeFraction);
            //pcout << "done" << std::endl << "Starting bubble history transition...";
            bubbleHistory->transition(*bubbleMatch, iT, bubbleVolumeRatio);
            //pcout << "done" << std::endl << "Starting bubble history update...";
            bubbleHistory->updateBubblePressure(fields.outsideDensity, param.rho_LB, param.alpha, param.beta);
            //pcout << "done" << std::endl;
        }

        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();
        global::timer("iteration").stop();

        if (param.useBubblePressureModel) {
            if (iT == 0 && param.freezeLargestBubble) {
                bubbleHistory->freezeLargestBubble();
            }
            if (iT % param.statIter == 0 || iT == param.maxIter-1) {
                bubbleHistory->timeHistoryLog(outDir + "bubbleTimeHistory.log");
                bubbleHistory->fullBubbleLog(outDir + "fullBubbleRecord.log");

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
                            iT * param.dt, it->first, v, r, center[0], center[1], center[2]);

                    // Save STL files of bubbles.

                    if (!util::fpequal(r, (T) 0.0, std::numeric_limits<T>::epsilon())) {
                        plint minNumOfTriangles = 64;
                        TriangleSet<T> triangleSet = constructSphere(center, r, minNumOfTriangles);
                        bubbles.append(triangleSet);
                    }
                }
                std::string fname = createFileName(outDir + "bubbles_", iT, PADDING) + ".stl";
                bubbles.writeBinarySTL(fname);
                pcout << "At iteration " << iT << ", the number of bubbles is " << numBubbles << std::endl;
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

    delete dynamics;

    exit(0);
}

