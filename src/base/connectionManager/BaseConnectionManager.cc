#include "BaseConnectionManager.h"

#include <cassert>

#include "NicEntryDebug.h"
#include "NicEntryDirect.h"
#include "BaseWorldUtility.h"
#include "FindModule.h"

#ifndef ccEV
#define ccEV (ev.isDisabled()||!coreDebug) ? ev : ev << getName() << ": "
#endif

namespace {
/**
 * On a torus the end and the begin of the axes are connected so you
 * get a circle. On a circle the distance between two points can't be greater
 * than half of the circumference.
 * If the normal distance between two points on one axis is bigger than
 * half of the size there must be a "shorter way" over the border on this axis
 */
static double dist(double coord1, double coord2, double size) {
    double difference = fabs(coord1 - coord2);
    if (difference == 0)
        // NOTE: event if size is zero
        return 0;
    else {
        assert(size != 0);
        double dist = FWMath::modulo(difference, size);
        return std::min(dist, size - dist);
    }
}

double sqrTorusDist(Coord c, Coord b, Coord size) {
    double xDist = dist(c.x, b.x, size.x);
    double yDist = dist(c.y, b.y, size.y);
    double zDist = dist(c.z, b.z, size.z);
    return xDist * xDist + yDist * yDist + zDist * zDist;
}
}

void BaseConnectionManager::initialize(int stage) {
    //BaseModule::initialize(stage);

    if (stage == 0) {
        if (hasPar("coreDebug"))
            coreDebug = par("coreDebug").boolValue();
        else
            coreDebug = false;

        drawMIR =
                hasPar("drawMaxIntfDist") ?
                        par("drawMaxIntfDist").boolValue() : false;

        ccEV << "initializing BaseConnectionManager\n";

        BaseWorldUtility* world =
                FindModule<BaseWorldUtility*>::findGlobalModule();

        assert(world != 0);

        playgroundSize = world->getPgs();
        useTorus = world->useTorus();

        if (hasPar("sendDirect"))
            sendDirect = par("sendDirect").boolValue();
        else
            sendDirect = false;

        // support VLC
        //maxInterferenceDistance = calcInterfDist();
        maxInterferenceDistance = 100;
        maxDistSquared = maxInterferenceDistance * maxInterferenceDistance;

        double maxHeadTXDistance = par("maxHeadTXDistance").doubleValue();
        double maxTailTXDistance = par("maxTailTXDistance").doubleValue();
        double maxHeadTXAngle = par("maxHeadTXAngle").doubleValue(); //degree
        double maxTailTXAngle = par("maxTailTXAngle").doubleValue(); //degree

        maxHeadTXAngle = M_PI*(maxHeadTXAngle/180); // rad
        maxTailTXAngle = M_PI*(maxTailTXAngle/180); // rad

        //----initialize node grid-----
        //step 1 - calculate dimension of grid
        //one cell should have at least the size of maxInterferenceDistance
        //but also should divide the playground in equal parts
        Coord dim((*playgroundSize) / maxInterferenceDistance);
        gridDim = GridCoord(dim);

        //A grid smaller or equal to 3x3 would mean that every cell has every
        //other cell as direct neighbor (if our playground is a torus, even if
        //not the most of the cells are direct neighbors of each other. So we
        //reduce the grid size to 1x1.
        if ((gridDim.x <= 3) && (gridDim.y <= 3) && (gridDim.z <= 3)) {
            gridDim.x = 1;
            gridDim.y = 1;
            gridDim.z = 1;
        } else {
            gridDim.x = std::max(1, gridDim.x);
            gridDim.y = std::max(1, gridDim.y);
            gridDim.z = std::max(1, gridDim.z);
        }

        //step 2 - initialize the matrix which represents our grid
        NicEntries entries;
        RowVector row;
        NicMatrix matrix;

        for (int i = 0; i < gridDim.z; ++i) {
            row.push_back(entries);         //copy empty NicEntries to RowVector
        }
        for (int i = 0; i < gridDim.y; ++i) { //fill the ColVector with copies of
            matrix.push_back(row);           //the RowVector.
        }
        for (int i = 0; i < gridDim.x; ++i) {   //fill the grid with copies of
            nicGrid.push_back(matrix);          //the matrix.
        }
        ccEV << " using " << gridDim.x << "x" << gridDim.y << "x" << gridDim.z
                    << " grid" << endl;

        //step 3 -  calculate the factor which maps the coordinate of a node
        //          to the grid cell
        //if we use a 1x1 grid every coordinate is mapped to (0,0, 0)
        findDistance = Coord(
                std::max(playgroundSize->x, maxInterferenceDistance),
                std::max(playgroundSize->y, maxInterferenceDistance),
                std::max(playgroundSize->z, maxInterferenceDistance));
        //otherwise we divide the playground into cells of size of the maximum
        //interference distance
        if (gridDim.x != 1)
            findDistance.x = playgroundSize->x / gridDim.x;
        if (gridDim.y != 1)
            findDistance.y = playgroundSize->y / gridDim.y;
        if (gridDim.z != 1)
            findDistance.z = playgroundSize->z / gridDim.z;

        //since the upper playground borders (at pg-size) are part of the
        //playground we have to assure that they are mapped to a valid
        //(the last) grid cell we do this by increasing the find distance
        //by a small value.
        //This also assures that findDistance is never zero.
        findDistance += Coord(EPSILON, EPSILON, EPSILON);

        //findDistance (equals cell size) has to be greater or equal
        //maxInt-distance
        assert(findDistance.x >= maxInterferenceDistance);
        assert(findDistance.y >= maxInterferenceDistance);
        assert(findDistance.z >= maxInterferenceDistance);

        //playGroundSize has to be part of the playGround
        assert(GridCoord(*playgroundSize, findDistance).x == gridDim.x - 1);
        assert(GridCoord(*playgroundSize, findDistance).y == gridDim.y - 1);
        assert(GridCoord(*playgroundSize, findDistance).z == gridDim.z - 1);
        ccEV << "findDistance is " << findDistance.info() << endl;
    } else if (stage == 1) {

    }
}

BaseConnectionManager::GridCoord BaseConnectionManager::getCellForCoordinate(
        const Coord& c) {
    return GridCoord(c, findDistance);
}

void BaseConnectionManager::updateConnections(int nicID, const Coord* oldPos,
        const Coord* newPos) {
    GridCoord oldCell = getCellForCoordinate(*oldPos);
    GridCoord newCell = getCellForCoordinate(*newPos);

    checkGrid(oldCell, newCell, nicID);
}

void BaseConnectionManager::updateConnectionsLite(int nicID,
        const Coord* oldPos, const Coord* newPos) {
    GridCoord oldCell = getCellForCoordinate(*oldPos);
    GridCoord newCell = getCellForCoordinate(*newPos);

    checkGridLite(oldCell, newCell, nicID);
}

BaseConnectionManager::NicEntries& BaseConnectionManager::getCellEntries(
        BaseConnectionManager::GridCoord& cell) {
    return nicGrid[cell.x][cell.y][cell.z];
}

void BaseConnectionManager::registerNicExt(int nicID) {
    NicEntries::mapped_type nicEntry = nics[nicID];

    GridCoord cell = getCellForCoordinate(nicEntry->pos);

    ccEV << " registering (ext) nic at loc " << cell.info() << std::endl;

    // add to matrix
    NicEntries& cellEntries = getCellEntries(cell);
    cellEntries[nicID] = nicEntry;
}

void BaseConnectionManager::checkGrid(BaseConnectionManager::GridCoord& oldCell,
        BaseConnectionManager::GridCoord& newCell, int id) {
    // structure to find union of grid squares
    CoordSet gridUnion(74);
    CoordSet gridUnionVLC(74);
    CoordSet gridUnionObstacle(74);

    // find nic at old position
    NicEntries& oldCellEntries = getCellEntries(oldCell);
    NicEntries::iterator it = oldCellEntries.find(id);
    NicEntries::mapped_type nic = it->second;

    // move nic to a new position in matrix
    if (oldCell != newCell) {
        oldCellEntries.erase(it);
        getCellEntries(newCell)[id] = nic;
    }

    if ((gridDim.x == 1) && (gridDim.y == 1) && (gridDim.z == 1)) {
        gridUnion.add(oldCell);
        gridUnionVLC.add(oldCell);
        gridUnionObstacle.add(oldCell);
    } else {
        //add grid around oldPos
        fillUnionWithNeighbors(gridUnion, oldCell);
        fillUnionWithNeighbors(gridUnionVLC, oldCell);
        fillUnionWithNeighbors(gridUnionObstacle, oldCell);

        if (oldCell != newCell) {
            //add grid around newPos
            fillUnionWithNeighbors(gridUnion, newCell);
            fillUnionWithNeighbors(gridUnionVLC, newCell);
            fillUnionWithNeighbors(gridUnionObstacle, newCell);
        }
    }

    ccEV << "updateNicConnections when position of nic #" << id << " updated."
                << endl;

    GridCoord* c = gridUnion.next();
    while (c != 0) {
        ccEV << "Update cons in [" << c->info() << "]" << endl;
        updateNicConnections(gridUnionVLC, gridUnionObstacle,
                getCellEntries(*c), nic);
        c = gridUnion.next();
    }
}

void BaseConnectionManager::checkGridLite(
        BaseConnectionManager::GridCoord& oldCell,
        BaseConnectionManager::GridCoord& newCell, int id) {
    // find nic at old position
    NicEntries& oldCellEntries = getCellEntries(oldCell);
    NicEntries::iterator it = oldCellEntries.find(id);
    NicEntries::mapped_type nic = it->second;

    // move nic to a new position in matrix
    if (oldCell != newCell) {
        oldCellEntries.erase(it);
        getCellEntries(newCell)[id] = nic;
    }
}

int BaseConnectionManager::wrapIfTorus(int value, int max) {
    if (value < 0) {
        if (useTorus) {
            return max + value;
        } else {
            return -1;
        }
    } else if (value >= max) {
        if (useTorus) {
            return value - max;
        } else {
            return -1;
        }
    } else {
        return value;
    }
}

void BaseConnectionManager::fillUnionWithNeighbors(CoordSet& gridUnion,
        GridCoord cell) {
    for (int iz = (int) cell.z - 1; iz <= (int) cell.z + 1; iz++) {
        int cz = wrapIfTorus(iz, gridDim.z);
        if (cz == -1) {
            continue;
        }
        for (int ix = (int) cell.x - 1; ix <= (int) cell.x + 1; ix++) {
            int cx = wrapIfTorus(ix, gridDim.x);
            if (cx == -1) {
                continue;
            }
            for (int iy = (int) cell.y - 1; iy <= (int) cell.y + 1; iy++) {
                int cy = wrapIfTorus(iy, gridDim.y);
                if (cy != -1) {
                    gridUnion.add(GridCoord(cx, cy, cz));
                }
            }
        }
    }
}

// is pToNic(rxNic) in range of pFromNic(txNic) ??
bool BaseConnectionManager::isInRange(CoordSet& gridUnionObstacle,
        BaseConnectionManager::NicEntries::mapped_type senderNic,
        BaseConnectionManager::NicEntries::mapped_type receiverNic) {

    int headlight = 1;
    int taillight = -1;

    ChannelAccess* senderModule = senderNic->chAccess;
    ChannelAccess* receiverModule = receiverNic->chAccess;

    // Position
    Coord senderPos = senderModule->getMobilityModule()->getCurrentPosition();
    Coord receiverPos =
            receiverModule->getMobilityModule()->getCurrentPosition();

    double distancefromSendertoReceiver = senderPos.distance(receiverPos);

    // Angle
    double senderAngle = senderModule->getMobilityModule()->getCurrentAngle();
    double receiverAngle =
            receiverModule->getMobilityModule()->getCurrentAngle();

    // Heading
    int senderHeading = senderModule->getMobilityModule()->getCurrentHeading();
    int receiverHeading =
            receiverModule->getMobilityModule()->getCurrentHeading();

    // Normalized Vector
    Coord vectorfromTXtoRX = (receiverPos - senderPos)
            / distancefromSendertoReceiver;
    Coord vectorTXheading = Coord(cos(senderAngle), sin(senderAngle))
            * senderHeading;
    Coord vectorRXheading = Coord(cos(receiverAngle), sin(receiverAngle))
            * receiverHeading;

    // skip the nic module on the same vehicle
    if (senderPos == receiverPos) {
        ccEV << "senderNic #" << senderNic->nicId << " and receiverNic #"
                    << receiverNic->nicId << " are on the same vehicle" << endl;
        return false;
    }

    // Sender is Headlight
    if (senderHeading == headlight) {
        ccEV << "senderNic #" << senderNic->nicId << " using Headlight" << endl;
        // Is it out of transmission distance?
        if (distancefromSendertoReceiver <= maxHeadTXDistance) {
            // Is it out of transmission angle?
            if ((vectorfromTXtoRX * vectorTXheading) > cos(maxHeadTXAngle)) {
                // Is it out of possible bearing?
                if ((vectorfromTXtoRX * vectorRXheading) < 0) {
                    // Is it blocked by vehicle(s)?
                    if (!isBlocked(gridUnionObstacle, senderNic, receiverNic)) {
                        ccEV << "TX/RX INFO : sender @ " << senderPos.info()
                                    << ", Angle = " << 180 * senderAngle / M_PI
                                    << ", Headlight?" << senderHeading
                                    << "; receiver @ " << receiverPos.info()
                                    << ", Angle = "
                                    << 180 * receiverAngle / M_PI
                                    << ", Headlight?" << receiverHeading
                                    << endl;

                        ccEV << "GOT RECEIVER : receverNic #"
                                    << receiverNic->nicId << " : Distance = "
                                    << distancefromSendertoReceiver
                                    << ", Angle = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI
                                    << ", Bearng to vectorfromTXtoRX = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << endl;
                        return true;
                    } else {
                        ccEV
                                    << "POSSIBLE RECEIVER BLOCKED BY VEHICLE(S) : receverNic #"
                                    << receiverNic->nicId << " : Distance = "
                                    << distancefromSendertoReceiver
                                    << ", Angle = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI
                                    << ", Bearng to vectorfromTXtoRX = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << endl;
                        return false;
                    }
                } else {
                    ccEV
                                << "OUT OF POSSIBLE BEARING to vectorfromTXtoRX (90, 180) or (-90, -180) ; receiverNic #"
                                << receiverNic->nicId
                                << " is on bearing to vectorfromTXtoRX -"
                                << 180
                                        * acos(
                                                vectorfromTXtoRX
                                                        * vectorRXheading)
                                        / M_PI << " or +"
                                << 180
                                        * acos(
                                                vectorfromTXtoRX
                                                        * vectorRXheading)
                                        / M_PI << " for senderNic #"
                                << senderNic->nicId << ", vectorfromTXtoRX ("
                                << vectorfromTXtoRX.x << ", "
                                << vectorfromTXtoRX.y << ") * vectorRXheading("
                                << vectorRXheading.x << ", "
                                << vectorRXheading.y << ") = "
                                << (vectorfromTXtoRX * vectorRXheading) << endl;
                    return false;
                }
            } else {
                ccEV
                            << "OUT OF POSSIBLE TRANSMISSION ANGLE (-45, 45) ; receiverNic #"
                            << receiverNic->nicId << " is on angle -"
                            << 180 * acos(vectorfromTXtoRX * vectorTXheading)
                                    / M_PI << " or +"
                            << 180 * acos(vectorfromTXtoRX * vectorTXheading)
                                    / M_PI << " of senderNic #"
                            << senderNic->nicId << ", vectorfromTXtoRX ("
                            << vectorfromTXtoRX.x << ", " << vectorfromTXtoRX.y
                            << ") * vectorTXheading(" << vectorTXheading.x
                            << ", " << vectorTXheading.y << ") = "
                            << (vectorfromTXtoRX * vectorTXheading) << endl;
                return false;
            }
        } else {
            ccEV
                        << "OUT OF POSSIBLE TRANSMISSION RANGE 100(m) ; distance between senderNic #"
                        << senderNic->nicId << " and receiverNic #"
                        << receiverNic->nicId << " = "
                        << distancefromSendertoReceiver << endl;
            return false;
        }
    }
    // Sender is Taillight
    else if (senderHeading == taillight) {
        ccEV << "senderNic #" << senderNic->nicId << " using Taillight" << endl;
        if (distancefromSendertoReceiver <= maxTailTXDistance) {
            if ((vectorfromTXtoRX * vectorTXheading) > cos(maxTailTXAngle)) {
                if ((vectorfromTXtoRX * vectorRXheading) < 0) {
                    if (!isBlocked(gridUnionObstacle, senderNic, receiverNic))
                        return true;
                    else {
                        ccEV
                                    << "POSSIBLE RECEIVER BLOCKED BY VEHICLE(S) : receverNic #"
                                    << receiverNic->nicId << " : Distance = "
                                    << distancefromSendertoReceiver
                                    << ", Angle = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorTXheading)
                                            / M_PI
                                    << ", Bearng to vectorfromTXtoRX = -"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << " or +"
                                    << 180
                                            * acos(
                                                    vectorfromTXtoRX
                                                            * vectorRXheading)
                                            / M_PI << endl;
                        return false;
                    }
                } else {
                    ccEV
                                << "OUT OF POSSIBLE BEARING to vectorfromTXtoRX (90, 180) or (-90, -180) ; receiverNic #"
                                << receiverNic->nicId
                                << " is on bearing to vectorfromTXtoRX -"
                                << 180
                                        * acos(
                                                vectorfromTXtoRX
                                                        * vectorRXheading)
                                        / M_PI << " or +"
                                << 180
                                        * acos(
                                                vectorfromTXtoRX
                                                        * vectorRXheading)
                                        / M_PI << " for senderNic #"
                                << senderNic->nicId << ", vectorfromTXtoRX ("
                                << vectorfromTXtoRX.x << ", "
                                << vectorfromTXtoRX.y << ") * vectorRXheading("
                                << vectorRXheading.x << ", "
                                << vectorRXheading.y << ") = "
                                << (vectorfromTXtoRX * vectorRXheading) << endl;
                    return false;
                }
            } else {
                ccEV
                            << "OUT OF POSSIBLE TRANSMISSION ANGLE (-60, 60) ; receiverNic #"
                            << receiverNic->nicId << " is on angle -"
                            << 180 * acos(vectorfromTXtoRX * vectorTXheading)
                                    / M_PI << " or +"
                            << 180 * acos(vectorfromTXtoRX * vectorTXheading)
                                    / M_PI << " of senderNic #"
                            << senderNic->nicId << ", vectorfromTXtoRX ("
                            << vectorfromTXtoRX.x << ", " << vectorfromTXtoRX.y
                            << ") * vectorTXheading(" << vectorTXheading.x
                            << ", " << vectorTXheading.y << ") = "
                            << (vectorfromTXtoRX * vectorTXheading) << endl;
                return false;
            }
        } else {
            ccEV
                        << "OUT OF POSSIBLE TRANSMISSION RANGE 30(m) ; distance between senderNic #"
                        << senderNic->nicId << " and receiverNic #"
                        << receiverNic->nicId << " = "
                        << distancefromSendertoReceiver << endl;
            return false;
        }
    } else
        ccEV << "unknown senderHeading" << senderHeading << endl;

    /*
     double dDistance = 0.0;

     if (useTorus) {
     dDistance = sqrTorusDist(pFromNic->pos, pToNic->pos, *playgroundSize);
     } else {
     dDistance = pFromNic->pos.sqrdist(pToNic->pos);
     }
     return (dDistance <= maxDistSquared);
     */
    return false;
}

bool BaseConnectionManager::isBlocked(CoordSet& gridUnionObstacle,
        BaseConnectionManager::NicEntries::mapped_type senderNic,
        BaseConnectionManager::NicEntries::mapped_type receiverNic) {

    GridCoord* cObstacle = gridUnionObstacle.next();
    while (cObstacle != 0) {
        NicEntries nmapObstacle = getCellEntries(*cObstacle);

        for (NicEntries::iterator k = nmapObstacle.begin();
                k != nmapObstacle.end(); ++k) {
            NicEntries::mapped_type obstacleNic = k->second;

            // no recursive connections
            if (senderNic->nicId == obstacleNic->nicId
                    || receiverNic->nicId == obstacleNic->nicId)
                continue;

            if (isObstacle(senderNic, receiverNic, obstacleNic)) {
                ccEV << "NLOS : Connection between senderNic #"
                            << senderNic->nicId << " and receiverNic #"
                            << receiverNic->nicId
                            << " is blocked by obstacleNic #"
                            << obstacleNic->nicId << endl;
                return true;
            }
        }

        cObstacle = gridUnionObstacle.next();
    }
    gridUnionObstacle.resetCurrent();

    return false;
}

bool BaseConnectionManager::isObstacle(
        BaseConnectionManager::NicEntries::mapped_type senderNic,
        BaseConnectionManager::NicEntries::mapped_type receiverNic,
        BaseConnectionManager::NicEntries::mapped_type obstacleNic) {

    double vehicleLength = 4.6; // 4620mm
    double vehicleWidth = 1.8; // 1775mm

    ChannelAccess* senderModule = senderNic->chAccess;
    ChannelAccess* receiverModule = receiverNic->chAccess;
    ChannelAccess* obstacleModule = obstacleNic->chAccess;

    // Position
    Coord senderPos = senderModule->getMobilityModule()->getCurrentPosition();
    Coord receiverPos =
            receiverModule->getMobilityModule()->getCurrentPosition();
    Coord obstaclePos =
            obstacleModule->getMobilityModule()->getCurrentPosition();

    // Angle
    double obstacleAngle =
            obstacleModule->getMobilityModule()->getCurrentAngle();

    // bounding box
    Coord boundBottomLeft = Coord(
            std::min(senderPos.x, receiverPos.x) - vehicleLength / 2,
            std::min(senderPos.y, receiverPos.y) - vehicleLength / 2);
    Coord boundUpperRight = Coord(
            std::max(senderPos.x, receiverPos.x) + vehicleLength / 2,
            std::max(senderPos.y, receiverPos.y) + vehicleLength / 2);

    if ((senderPos == obstaclePos) || (receiverPos == obstaclePos)) {
        ccEV << "Skip possible obstacleNic #" << obstacleNic->nicId
                    << " is installed on either senderNic #" << senderNic->nicId
                    << "or receiverNic #" << receiverNic->nicId << "'s vehicle"
                    << endl;
        return false;
    }

    // vehicles located at point where out of bounding box has no chance to block connection between sender and receiver
    if (obstaclePos.x < boundBottomLeft.x || obstaclePos.y < boundBottomLeft.y
            || obstaclePos.x > boundUpperRight.x
            || obstaclePos.y > boundUpperRight.y)
        return false;
    else {
        // y=mx+b
        double m = (receiverPos.y - senderPos.y)
                / (receiverPos.x - senderPos.x);
        double b = receiverPos.y - m * receiverPos.x;

        // four corner position of vehicle
        Coord vehicleCorner1 = obstaclePos
                + Coord(
                        cos(obstacleAngle) * (vehicleLength / 2)
                                - sin(obstacleAngle) * (vehicleWidth / 2),
                        sin(obstacleAngle) * (vehicleLength / 2)
                                + cos(obstacleAngle) * (vehicleWidth / 2));
        Coord vehicleCorner2 = obstaclePos
                + Coord(
                        cos(obstacleAngle) * (-vehicleLength / 2)
                                - sin(obstacleAngle) * (vehicleWidth / 2),
                        sin(obstacleAngle) * (-vehicleLength / 2)
                                + cos(obstacleAngle) * (vehicleWidth / 2));
        Coord vehicleCorner3 = obstaclePos
                + Coord(
                        cos(obstacleAngle) * (-vehicleLength / 2)
                                - sin(obstacleAngle) * (-vehicleWidth / 2),
                        sin(obstacleAngle) * (-vehicleLength / 2)
                                + cos(obstacleAngle) * (-vehicleWidth / 2));
        Coord vehicleCorner4 = obstaclePos
                + Coord(
                        cos(obstacleAngle) * (vehicleLength / 2)
                                - sin(obstacleAngle) * (-vehicleWidth / 2),
                        sin(obstacleAngle) * (vehicleLength / 2)
                                + cos(obstacleAngle) * (-vehicleWidth / 2));

        // f(x,y)= mx-y+b > 0 Left-side of line ; < 0 right-side
        if (((m * vehicleCorner1.x - vehicleCorner1.y + b) > 0)
                && ((m * vehicleCorner2.x - vehicleCorner2.y + b) > 0)
                && ((m * vehicleCorner3.x - vehicleCorner3.y + b) > 0)
                && ((m * vehicleCorner4.x - vehicleCorner4.y + b) > 0)) {
            // left-side of f(x,y)
            return false;
        } else if (((m * vehicleCorner1.x - vehicleCorner1.y + b) < 0)
                && ((m * vehicleCorner2.x - vehicleCorner2.y + b) < 0)
                && ((m * vehicleCorner3.x - vehicleCorner3.y + b) < 0)
                && ((m * vehicleCorner4.x - vehicleCorner4.y + b) < 0)) {
            // right-side of f(x,y)
            return false;
        } else {
            ccEV << "senderPos = " << senderPos.info() << ", receiverPos = "
                        << receiverPos.info() << endl;

            ccEV << "boundBottomLeft = " << boundBottomLeft.info()
                        << ", boundUpperRight = " << boundUpperRight.info()
                        << endl;

            ccEV << "y = " << m << "x + (" << b << ") ; obstacleNic #"
                        << obstacleNic->nicId << "whose position = "
                        << obstaclePos.info() << ", angle = "
                        << 180 * obstacleAngle / M_PI << "(degree) = "
                        << obstacleAngle << "(rad) and cornersPos = "
                        << vehicleCorner1.info() << " " << vehicleCorner2.info()
                        << " " << vehicleCorner3.info() << " "
                        << vehicleCorner4.info() << endl;

            ccEV << "cos(angle) = " << cos(obstacleAngle) << ", sin(angle) = "
                        << sin(obstacleAngle) << endl;

            return true;
        }
    }
}

// Support filter out NLOS blocked by vehicles
/*
 void BaseConnectionManager::updateNicConnections(NicEntries& nmap,
 BaseConnectionManager::NicEntries::mapped_type nic) {
 int id = nic->nicId;

 for (NicEntries::iterator i = nmap.begin(); i != nmap.end(); ++i) {
 NicEntries::mapped_type nic_i = i->second;

 // no recursive connections
 if (nic_i->nicId == id)
 continue;

 bool inRange = isInRange(nic, nic_i);
 bool connected = nic->isConnected(nic_i);

 if (inRange && !connected) {
 // nodes within communication range: connect
 // nodes within communication range && not yet connected
 ccEV << "nic #" << id << " and #" << nic_i->nicId << " are in range"
 << endl;
 nic->connectTo(nic_i);
 nic_i->connectTo(nic);
 } else if (!inRange && connected) {
 // out of range: disconnect
 // out of range, and still connected
 ccEV << "nic #" << id << " and #" << nic_i->nicId
 << " are NOT in range" << endl;
 nic->disconnectFrom(nic_i);
 nic_i->disconnectFrom(nic);
 }
 }
 }
 */

void BaseConnectionManager::updateNicConnections(CoordSet& gridUnionVLC,
        CoordSet& gridUnionObstacle, NicEntries& nmap,
        BaseConnectionManager::NicEntries::mapped_type nic) {

    for (NicEntries::iterator i = nmap.begin(); i != nmap.end(); ++i) {
        NicEntries::mapped_type nic_tx = i->second;

        GridCoord* cVLC = gridUnionVLC.next();
        while (cVLC != 0) {
            ccEV << "vs cons in [" << cVLC->info() << "]" << endl;
            NicEntries nmapVLC = getCellEntries(*cVLC);

            for (NicEntries::iterator j = nmapVLC.begin(); j != nmapVLC.end();
                    ++j) {
                NicEntries::mapped_type nic_rx = j->second;

                // no recursive connections
                if (nic_tx->nicId == nic_rx->nicId)
                    continue;

                // is nic_rx in range of nic_tx ?
                bool inRange = isInRange(gridUnionObstacle, nic_tx, nic_rx);
                bool connected = nic_tx->isConnected(nic_rx);

                if (inRange && !connected) {
                    // nodes within communication range: connect
                    // nodes within communication range && not yet connected
                    ccEV << "nic_rx #" << nic_rx->nicId
                                << " are in range of nic_tx #" << nic_tx->nicId
                                << endl;
                    nic_tx->connectTo(nic_rx);
                } else if (!inRange && connected) {
                    // out of range: disconnect
                    // out of range, and still connected
                    ccEV << "nic_rx #" << nic_rx->nicId
                                << " are NOT in range of nic_tx # "
                                << nic_tx->nicId << endl;
                    nic_tx->disconnectFrom(nic_rx);
                }
            }
            cVLC = gridUnionVLC.next();
        }
        gridUnionVLC.resetCurrent();
    }
}

bool BaseConnectionManager::registerNic(cModule* nic, ChannelAccess* chAccess,
        const Coord* nicPos) {
    assert(nic != 0);

    int nicID = nic->getId();
    ccEV << " registering nic #" << nicID << endl;

// create new NicEntry
    NicEntries::mapped_type nicEntry;

    if (sendDirect)
        nicEntry = new NicEntryDirect(coreDebug);
    else
        nicEntry = new NicEntryDebug(coreDebug);

// fill nicEntry
    nicEntry->nicPtr = nic;
    nicEntry->nicId = nicID;
    nicEntry->hostId = nic->getParentModule()->getId();
    nicEntry->pos = *nicPos;
    nicEntry->chAccess = chAccess;

// add to map
    nics[nicID] = nicEntry;

    registerNicExt(nicID);

    updateConnections(nicID, nicPos, nicPos);

    if (drawMIR) {
        nic->getParentModule()->getDisplayString().setTagArg("r", 0,
                maxInterferenceDistance);
    }

    return sendDirect;
}

bool BaseConnectionManager::registerNicLite(cModule* nic,
        ChannelAccess* chAccess, const Coord* nicPos) {
    assert(nic != 0);

    int nicID = nic->getId();
    ccEV << " registering nic #" << nicID << endl;

// create new NicEntry
    NicEntries::mapped_type nicEntry;

    if (sendDirect)
        nicEntry = new NicEntryDirect(coreDebug);
    else
        nicEntry = new NicEntryDebug(coreDebug);

// fill nicEntry
    nicEntry->nicPtr = nic;
    nicEntry->nicId = nicID;
    nicEntry->hostId = nic->getParentModule()->getId();
    nicEntry->pos = *nicPos;
    nicEntry->chAccess = chAccess;

// add to map
    nics[nicID] = nicEntry;

    registerNicExt(nicID);

    updateConnections(nicID, nicPos, nicPos);

    if (drawMIR) {
        nic->getParentModule()->getDisplayString().setTagArg("r", 0,
                maxInterferenceDistance);
    }

    return sendDirect;
}

bool BaseConnectionManager::unregisterNic(cModule* nicModule) {
    assert(nicModule != 0);

// find nicEntry
    int nicID = nicModule->getId();
    ccEV << " unregistering nic #" << nicID << endl;

//we assume that the module was previously registered with this CM
//TODO: maybe change this to an omnet-error instead of an assertion
    assert(nics.find(nicID) != nics.end());
    NicEntries::mapped_type nicEntry = nics[nicID];

// get all affected grid squares
    CoordSet gridUnion(74);
    GridCoord cell = getCellForCoordinate(nicEntry->pos);
    if ((gridDim.x == 1) && (gridDim.y == 1) && (gridDim.z == 1)) {
        gridUnion.add(cell);
    } else {
        fillUnionWithNeighbors(gridUnion, cell);
    }

// disconnect from all NICs in these grid squares
    GridCoord* c = gridUnion.next();
    while (c != 0) {
        ccEV << "Update cons in [" << c->info() << "]" << endl;
        NicEntries& nmap = getCellEntries(*c);
        for (NicEntries::iterator i = nmap.begin(); i != nmap.end(); ++i) {
            NicEntries::mapped_type other = i->second;
            if (other == nicEntry)
                continue;

            if (nicEntry->isConnected(other))
                nicEntry->disconnectFrom(other);

            if (other->isConnected(nicEntry))
                other->disconnectFrom(nicEntry);
        }
        c = gridUnion.next();
    }

// erase from grid
    NicEntries& cellEntries = getCellEntries(cell);
    cellEntries.erase(nicID);

// erase from list of known nics
    nics.erase(nicID);

    delete nicEntry;

    return true;
}

void BaseConnectionManager::updateNicPos(int nicID, const Coord* newPos) {
    NicEntries::iterator ItNic = nics.find(nicID);
    if (ItNic == nics.end())
        error(
                "No nic with this ID (%d) is registered with this ConnectionManager.",
                nicID);

    Coord oldPos = ItNic->second->pos;
    ItNic->second->pos = *newPos;

    updateConnections(nicID, &oldPos, newPos);
}

void BaseConnectionManager::updateNicPosLite(int nicID, const Coord* newPos) {
    NicEntries::iterator ItNic = nics.find(nicID);
    if (ItNic == nics.end())
        error(
                "No nic with this ID (%d) is registered with this ConnectionManager.",
                nicID);

    Coord oldPos = ItNic->second->pos;
    ItNic->second->pos = *newPos;

    updateConnectionsLite(nicID, &oldPos, newPos);
}

const NicEntry::GateList& BaseConnectionManager::getGateList(int nicID) const {
    NicEntries::const_iterator ItNic = nics.find(nicID);
    if (ItNic == nics.end())
        error(
                "No nic with this ID (%d) is registered with this ConnectionManager.",
                nicID);

    return ItNic->second->getGateList();
}

const cGate* BaseConnectionManager::getOutGateTo(const NicEntry* nic,
        const NicEntry* targetNic) const {
    NicEntries::const_iterator ItNic = nics.find(nic->nicId);
    if (ItNic == nics.end())
        error(
                "No nic with this ID (%d) is registered with this ConnectionManager.",
                nic->nicId);

    return ItNic->second->getOutGateTo(targetNic);
}

BaseConnectionManager::~BaseConnectionManager() {
    for (NicEntries::iterator ne = nics.begin(); ne != nics.end(); ne++) {
        delete ne->second;
    }
}
