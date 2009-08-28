/*
 *  mwMaskStimulusFactory.h
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */


#ifndef mwMaskStimulus_FACTORY_H
#define mwMaskStimulus_FACTORY_H

#include "mwMaskStimulus.h"

#include "MonkeyWorksCore/ComponentFactory.h"
using namespace mw;

class mwMaskStimulusFactory : public ComponentFactory {
	virtual shared_ptr<mw::Component> createObject(std::map<std::string, std::string> parameters,
                                                   mw::ComponentRegistry *reg);
};

#endif
