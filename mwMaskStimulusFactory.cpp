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
    
    // generate default seed
    // BJG: this is sort of a bad way to generate the seed as it only changes every second
    uint32_t seed = static_cast<unsigned int>(std::time(0));
    std::stringstream out;
    std::string seedStr;
    out << seed;
    seedStr = out.str();
    
    // get random seed from 
    //shared_ptr<Variable> random_seed = reg->getVariable(parameters["random_seed"], seedStr);
    // TODO: set this to some default but unique value, like the current time in microseconds
    uint32_t random_seed = seed;
    if(!parameters["random_seed"].empty()){
        random_seed = (uint32_t)(reg->getNumber(parameters["random_seed"]).getFloat());
        // TODO error checking
        mprintf("Found random_seed of %i",random_seed);
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
                                                                                              random_seed));
	//shared_ptr <mwMaskStimulus> newMaskStimulus = shared_ptr<mwMaskStimulus>(new mwMaskStimulus(tagname, another_attribute));

    bool deferred = true;
    // bjg: do not read the deferred varialbe from XML, it should always be true
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
