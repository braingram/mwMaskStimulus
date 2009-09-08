/*
 *  mwMaskStimulusFactory.cpp
 *  mwMaskStimulus
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */

#include "mwMaskStimulusFactory.h"

#include <boost/regex.hpp>
#include <boost/filesystem.hpp>
#include <sys/time.h>
#include <sstream>

#include "MonkeyWorksCore/ComponentRegistry.h"

using namespace mw;

shared_ptr<mw::Component> mwMaskStimulusFactory::createObject(std::map<std::string, std::string> parameters,
                                                               mw::ComponentRegistry *reg) {
	
    // ----------------------
    // --- Get Attributes ---
    // ----------------------
    
    // alpha is not required, and uses a default value below
    // same with seed (to type uint32_t)
    REQUIRE_ATTRIBUTES(parameters, 
					   "tag",
                       "x_size",
                       "y_size",
                       "x_position",
                       "y_position",
                       "rotation",
                       "path");

    // standard 'tag' of the component, think of this as the 'name'
	std::string tagname(parameters.find("tag")->second);

    // other variables of the component go here
    // for a mask I want:
    //   current random seed
    //   random program used (boost::random, merzian(sp?))
    //   (next few taken from stimulus)
    //   x_size, y_size
    //   x_position, y_position
    //   rotation
    //   path
	shared_ptr<Variable> x_size = reg->getVariable(parameters.find("x_size")->second);
	shared_ptr<Variable> y_size = reg->getVariable(parameters.find("y_size")->second);
	shared_ptr<Variable> x_position = reg->getVariable(parameters.find("x_position")->second);
	shared_ptr<Variable> y_position = reg->getVariable(parameters.find("y_position")->second);
	shared_ptr<Variable> rotation = reg->getVariable(parameters.find("rotation")->second);
    // alpha_multiplier has a default value of 1.0
    shared_ptr<Variable> alpha_multiplier = reg->getVariable(parameters["alpha_multiplier"], std::string("1.0"));
    
//    // generate default seed
//    shared_ptr<Variable> random_seed = reg->getVariable(parameters.find("random_seed")->second);
//    std::cout << random_seed << "\n";
//    if (random_seed == 0) {
//        // no random_seed was defined, so generate one
//        // BJG: this is sort of a bad way to generate the seed as it only changes every second 
//        // TODO: set this to some default but unique value, like the current time in microseconds
//        uint32_t default_seed = static_cast<unsigned int>(std::time(0));
//        random_seed->setValue(Data((long)default_seed));
//    }
    // get random seed from 
    //shared_ptr<Variable> random_seed = reg->getVariable(parameters["random_seed"], seedStr);
    struct timeval tv;
    gettimeofday(&tv, NULL);
    uint32_t default_seed = (tv.tv_sec % 1000) * 1000000 + tv.tv_usec;
    
    //uint32_t default_seed = static_cast<unsigned int>(std::time(0));
    shared_ptr<Variable> random_seed(new ConstantVariable(Data((long)default_seed))); 
    if(!parameters["random_seed"].empty()){
        random_seed = reg->getVariable(parameters.find("random_seed")->second);
        //random_seed = (uint32_t)(reg->getNumber(parameters["random_seed"]).getFloat());
        // TODO error checking
        mprintf("Found random_seed of %i",random_seed->getValue().getInteger());
    }
    
//    shared_ptr<Variable> random_phase_per_channel = reg->getVariable(parameters.find("random_phase_per_channel")->second);
//    if (random_phase_per_channel == 0) {
//        random_phase_per_channel->setValue(Data(false));
//    }
    shared_ptr<Variable> random_phase_per_channel(new ConstantVariable(Data(false)));
    if (!parameters["random_phase_per_channel"].empty()) {
        random_phase_per_channel = reg->getVariable(parameters.find("random_seed")->second);
    }
    //long random_seed = reg->getLong(parameters["random_seed"], seedStr);
	// !!! check this next one !!!
    boost::filesystem::path full_path = reg->getPath(parameters["working_path"], parameters["path"]);
    
    
    // ----------------------
    // --- Error Checking ---
    // ----------------------
    
	checkAttribute(x_size, parameters.find("reference_id")->second, "x_size", parameters["x_size"]);
	checkAttribute(y_size, parameters.find("reference_id")->second, "y_size", parameters.find("y_size")->second);
	checkAttribute(x_position, parameters.find("reference_id")->second, "x_position", parameters.find("x_position")->second);
	checkAttribute(y_position, parameters.find("reference_id")->second, "y_position", parameters.find("y_position")->second);
	checkAttribute(rotation, parameters.find("reference_id")->second, "rotation", parameters.find("rotation")->second);
	checkAttribute(alpha_multiplier,parameters.find("reference_id")->second, "alpha_multiplier", parameters.find("alpha_multiplier")->second);
    checkAttribute(random_seed, parameters.find("reference_id")->second, "random_seed", parameters.find("random_seed")->second);
    checkAttribute(random_phase_per_channel, parameters.find("reference_id")->second, "random_phase_per_channel", parameters.find("random_phase_per_channel")->second);
	if(!boost::filesystem::exists(full_path)) {
		throw InvalidReferenceException(parameters.find("reference_id")->second, "path", parameters.find("path")->second);
	}
	
	if(boost::filesystem::is_directory(full_path)) {
		throw InvalidReferenceException(parameters.find("reference_id")->second, "path", parameters.find("path")->second);
	}
	
	if(GlobalCurrentExperiment == NULL) {
		throw SimpleException("no experiment currently defined");		
	}
	
	shared_ptr<StimulusDisplay> defaultDisplay = GlobalCurrentExperiment->getStimulusDisplay();
	if(defaultDisplay == 0) {
		throw SimpleException("no stimulusDisplay in current experiment");
	}
    
	//shared_ptr<Variable> another_attribute = reg->getVariable(parameters.find("another_attribute")->second);
	

    // ----------------------
    // --- Make Component ---
    // ----------------------
    
    shared_ptr <mwMaskStimulus> newMaskStimulus = shared_ptr<mwMaskStimulus>(new mwMaskStimulus(tagname,
                                                                                              full_path.string(),
                                                                                              x_position,
                                                                                              y_position,
                                                                                              x_size,
                                                                                              y_size,
                                                                                              rotation,
                                                                                              alpha_multiplier,
                                                                                              random_seed,
                                                                                              random_phase_per_channel));
	//shared_ptr <mwMaskStimulus> newMaskStimulus = shared_ptr<mwMaskStimulus>(new mwMaskStimulus(tagname, another_attribute));
    
    // bjg: do not read the deferred varialbe from XML, it should always be true
    bool deferred = true;
    //if(!parameters["deferred"].empty()){
    //    deferred = reg->getBoolean(parameters["deferred"]);
    //}
    
    // TODO: deferred load?
    if(!deferred){
        newMaskStimulus->load(defaultDisplay.get());
    }
    
	shared_ptr <StimulusNode> thisStimNode = shared_ptr<StimulusNode>(new StimulusNode(newMaskStimulus));
	reg->registerStimulusNode(tagname, thisStimNode);
	
	//return newImageStimulus;
	return newMaskStimulus;
}
