/***************************************************************************
 * file:        ChannelAccess.cc
 *
 * author:      Marc Loebbers
 *

 * copyright:   (C) 2004 Telecommunication Networks Group (TKN) at
 *              Technische Universitaet Berlin, Germany.
 *
 *              This program is free software; you can redistribute it
 *              and/or modify it under the terms of the GNU General Public
 *              License as published by the Free Software Foundation; either
 *              version 2 of the License, or (at your option) any later
 *              version.
 *              For further information see file COPYING
 *              in the top level directory
 ***************************************************************************
 * part of:     framework implementation developed by tkn
 * description: - Base class for physical layers
 *              - if you create your own physical layer, please subclass
 *                from this class and use the sendToChannel() function!!
 ***************************************************************************
 * changelog:   $Revision: 284 $
 *              last modified:   $Date: 2006-06-07 16:55:24 +0200 (Mi, 07 Jun 2006) $
 *              by:              $Author: willkomm $
 **************************************************************************/

#include "base/connectionManager/ChannelAccess.h"

#include <cassert>

#include "FindModule.h"
#include "BaseWorldUtility.h"
#include "BaseConnectionManager.h"

// for Hello Message
#include "NicEntryDebug.h"
#include "NicEntryDirect.h"
#include <AirFrame_m.h>

using std::endl;

const simsignalwrap_t ChannelAccess::mobilityStateChangedSignal =
        simsignalwrap_t(MIXIM_SIGNAL_MOBILITY_CHANGE_NAME);

BaseConnectionManager* ChannelAccess::getConnectionManager(cModule* nic) {
    std::string cmName =
            nic->hasPar("connectionManagerName") ?
                    nic->par("connectionManagerName").stringValue() : "";
    if (cmName != "") {
        cModule* ccModule = simulation.getModuleByPath(cmName.c_str());

        return dynamic_cast<BaseConnectionManager *>(ccModule);
    } else {
        return FindModule<BaseConnectionManager *>::findGlobalModule();
    }
}

void ChannelAccess::initialize(int stage) {
    BatteryAccess::initialize(stage);

    if (stage == 0) {
        hasPar("coreDebug") ?
                coreDebug = par("coreDebug").boolValue() : coreDebug = false;

        findHost()->subscribe(mobilityStateChangedSignal, this);

        cModule* nic = getParentModule();
        cc = getConnectionManager(nic);
        if (cc == NULL)
            error("Could not find connectionmanager module");
        isRegistered = false;
    }

    usePropagationDelay = par("usePropagationDelay");
}

void ChannelAccess::sendToChannel(cPacket *msg) {
    const NicEntry::GateList& gateList = cc->getGateList(
            getParentModule()->getId());
    NicEntry::GateList::const_iterator i = gateList.begin();

    AirFrame* airframe = static_cast<AirFrame*>(msg);
    bool asymmetry;
    if (useSendDirect) {
        // use Andras stuff
        if (i != gateList.end()) {
            simtime_t delay = SIMTIME_ZERO;
            for (; i != --gateList.end(); ++i) {

                // nic->rx = i->first : nic_tx -> nic_rx must exist
                // and isAsymmetry(nic_txID, nic_rxID) tell us nic->rx -> nic_tx exist or not
                // if not, isAsymmetry will return true
                asymmetry = cc->isAsymmetry(getParentModule()->getId(),
                        i->first->nicId);

                if (asymmetry) {
                    coreEV
                                  << "nic_rx doesn't connect to nic_tx --> Asymmetry happen."
                                  << endl;
                    airframe->setIsAsymmetry(true);
                } else {
                    coreEV << "nic_tx <-> nic_rx --> Symmetry" << endl;
                    airframe->setIsAsymmetry(false);
                }

                //calculate delay (Propagation) to this receiving nic
                delay = calculatePropagationDelay(i->first);

                int radioStart = i->second->getId();
                int radioEnd = radioStart + i->second->size();
                for (int g = radioStart; g != radioEnd; ++g)
                    sendDirect(static_cast<cPacket*>(msg->dup()), delay,
                            msg->getDuration(), i->second->getOwnerModule(), g);
            }

            asymmetry = cc->isAsymmetry(getParentModule()->getId(),
                    i->first->nicId);

            if (asymmetry) {
                coreEV
                              << "nic_rx doesn't connect to nic_tx --> Asymmetry happen."
                              << endl;
                airframe->setIsAsymmetry(true);
            } else {
                coreEV << "nic_tx <-> nic_rx --> Symmetry" << endl;
                airframe->setIsAsymmetry(false);
            }

            //calculate delay (Propagation) to this receiving nic
            delay = calculatePropagationDelay(i->first);

            int radioStart = i->second->getId();
            int radioEnd = radioStart + i->second->size();
            for (int g = radioStart; g != --radioEnd; ++g)
                sendDirect(static_cast<cPacket*>(msg->dup()), delay,
                        msg->getDuration(), i->second->getOwnerModule(), g);

            sendDirect(msg, delay, msg->getDuration(),
                    i->second->getOwnerModule(), radioEnd);

        } else {
            coreEV << "Nic is not connected to any gates!" << endl;
            delete msg;
        }
    } else {
        // use our stuff
        coreEV << "sendToChannel: sending to gates\n";
        if (i != gateList.end()) {
            simtime_t delay = SIMTIME_ZERO;
            for (; i != --gateList.end(); ++i) {
                //calculate delay (Propagation) to this receiving nic
                delay = calculatePropagationDelay(i->first);

                sendDelayed(static_cast<cPacket*>(msg->dup()), delay,
                        i->second);
            }
            //calculate delay (Propagation) to this receiving nic
            delay = calculatePropagationDelay(i->first);

            sendDelayed(msg, delay, i->second);
        } else {
            coreEV << "Nic is not connected to any gates!" << endl;
            delete msg;
        }
    }
}

simtime_t ChannelAccess::calculatePropagationDelay(const NicEntry* nic) {
    if (!usePropagationDelay)
        return 0;

    ChannelAccess * const senderModule = this;
    ChannelAccess * const receiverModule = nic->chAccess;
    //const simtime_t_cref sStart         = simTime();

    assert(senderModule);
    assert(receiverModule);

    /** claim the Move pattern of the sender from the Signal */
    Coord sendersPos =
            senderModule->getMobilityModule()->getCurrentPosition(/*sStart*/);
    Coord receiverPos =
            receiverModule->getMobilityModule()->getCurrentPosition(/*sStart*/);

    // this time-point is used to calculate the distance between sending and receiving host
    return receiverPos.distance(sendersPos) / BaseWorldUtility::speedOfLight;
}

void ChannelAccess::receiveSignal(cComponent *source, simsignal_t signalID,
        cObject *obj) {
    if (signalID == mobilityStateChangedSignal) {
        ChannelMobilityPtrType const mobility = check_and_cast<
                ChannelMobilityPtrType>(obj);
        Coord pos = mobility->getCurrentPosition();
        /*
         if (isRegistered) {
         cc->updateNicPos(getParentModule()->getId(), &pos);
         } else {
         // register the nic with ConnectionManager
         // returns true, if sendDirect is used
         useSendDirect = cc->registerNic(getParentModule(), this, &pos);
         isRegistered = true;
         }
         */

        if (this->getMobilityModule()->getCurrentHeading() == -1) {
            coreEV << "mobilityStateChangedSignal for tail module" << endl;
            if (isRegistered) {
                cc->updateNicPos(getParentModule()->getId(), &pos);
            } else {
                // register the nic with ConnectionManager
                // returns true, if sendDirect is used
                useSendDirect = cc->registerNic(getParentModule(), this, &pos);
                isRegistered = true;
            }
        } else {
            coreEV << "mobilityStateChangedSignal for head module" << endl;
            if (isRegistered) {
                cc->updateNicPosLite(getParentModule()->getId(), &pos);
            } else {
                // register the nic with ConnectionManager
                // returns true, if sendDirect is used
                useSendDirect = cc->registerNicLite(getParentModule(), this,
                        &pos);
                isRegistered = true;
            }
        }
    }
}
