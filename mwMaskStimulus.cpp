/*
 *  mwMaskStimulus.cpp
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */

#include "mwMaskStimulus.h"


mwMaskStimulus::mwMaskStimulus(std::string _tag, std::string _filename,
                               shared_ptr<Variable> _xoffset,
                               shared_ptr<Variable> _yoffset,
                               shared_ptr<Variable> _xscale,
                               shared_ptr<Variable> _yscale,
                               shared_ptr<Variable> _rot,
                               shared_ptr<Variable> _alpha,
                               shared_ptr<Variable> _random_seed)
                                : ImageStimulus (_tag, _filename, _xoffset, _yoffset, _xscale, _yscale, _rot, _alpha) {
    //mask_textures = new ExpandableList<GLuint>();
    random_seed = _random_seed;
}

mwMaskStimulus::mwMaskStimulus(const mwMaskStimulus &tocopy)
    :ImageStimulus((ImageStimulus&) tocopy) {
        random_seed = tocopy.random_seed;
}

mwMaskStimulus::~mwMaskStimulus(){
    if(loaded){
        free(image_data);
    }
}

shared_ptr<Variable> mwMaskStimulus::getRandomSeed() {
    return random_seed;
}

void mwMaskStimulus::makeMask(StimulusDisplay *display) {
    if(!loaded){
        return;
    }
    
    // do fourier transform
      
    
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//    // copy new mask to textures
//    for(int i = 0; i < display->getNContexts(); i++){
//        display->setCurrent(i);
//        
//    };
}

void mwMaskStimulus::load(StimulusDisplay *display) {
    // loads image into a texture of format RGBA FLOAT
    
    //GlobalOpenGLContextManager->setGlobalContextCurrent();
    if(loaded){
		return;
	}
	
    // ---------- New ---------------
    
    //std::cout << "Loading image mask stimulus\n";
    // check if file exists
    if(filename == ""){
		throw SimpleException("Cannot load image (NULL filename).");
		return;
	}
    FILE *test = fopen(filename.c_str(),"r");
	if(!test){
		throw SimpleException("Image file does not exist", filename);
	}
	fclose(test);
    
    // start up Devil
    ilInit();
    ilutInit(); // !!! necessary?
    ilutRenderer(ILUT_OPENGL); // !!! necessary?
    ilutEnable(ILUT_OPENGL_CONV); // !!! necessary?
    
    //load image from FILE
    GLuint texture_map;
    ILuint il_image_name;
	ILenum il_error;
    ilGenImages(1,&il_image_name);
    ilBindImage(il_image_name);
    ilLoadImage(filename.c_str());
    ilConvertImage(IL_RGBA, IL_FLOAT);
    if((il_error = ilGetError()) != IL_NO_ERROR) {
        // TODO HANDLE ERROR
        merror(M_DISPLAY_MESSAGE_DOMAIN,
               "IL Image Library Error: %x", 
               il_error);
		
		throw SimpleException("Cannot load image", filename);
		//lock->unlock();
        return;
    }
    
    width = (int)ilGetInteger(IL_IMAGE_WIDTH);
    height = (int)ilGetInteger(IL_IMAGE_HEIGHT);
    image_data = (float*) calloc(4 * height * width, sizeof(float));
    
    // copy data to image_data
    ilCopyPixels(0, 0, 0, (ILuint)width, (ILuint)height, 1, IL_RGBA, IL_FLOAT, (ILvoid*) image_data);
    // !!! TODO check for errors
    if ((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(
                              (boost::format("Cannot copy image data from ilImage") % il_error).str());
        //image_loaded = false;
        free(image_data);
        return;
    }
    //ILubyte *Data = ilGetData(); 
    
    // move 'masks' (original images right now) to gpu
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    for(int i = 0; i < display->getNContexts(); i++){
		display->setCurrent(i);        
		texture_map = ilutGLBindMipmaps();
        texture_maps->addElement(i, texture_map);
		if(texture_map){
			mprintf("Image loaded into texture_map %d", texture_map);
		}
	}
    
    if((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(
                              (boost::format("Cannot bind image texture. IL Image Library Error") % il_error).str());
        //image_loaded = false;
		//lock->unlock();
        free(image_data);
        return;
    }
	
	ilDeleteImage(il_image_name);
	if((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(M_DISPLAY_MESSAGE_DOMAIN,
                              (boost::format("Cannot delete image. IL Image Library Error: %s") % il_error).str());
        //image_loaded = false;
		//lock->unlock();
        free(image_data);
        return;
    }
    
    loaded = true;
    
    // make mask
    //makeMask(display);
    
    // -------- end New -------------
    
//	// TODO: this needs clean up.  We are counting on all of the contexts
//	// in the stimulus display to have the exact same history.  Ideally, this
//	// should be true, but we should eventually be robust in case it isn't
//	for(int i = 0; i < display->getNContexts(); i++){
//		display->setCurrent(i);
//		GLuint texture_map = OpenGLImageLoader::load(filename, display, &width, &height);
//		texture_maps->addElement(i, texture_map);
//        //		fprintf(stderr, "Loaded texture map %u into context %d\n", (unsigned int)texture_map, i);fflush(stderr);
//		if(texture_map){
//			mprintf("Image loaded into texture_map %d", texture_map);
//		}
//	}
//    
//    // BJG: is it necessary to store CPU-side image data for each context?
//    if (display->getNContexts() < 1) {
//        throw SimpleException("no display contexts currently defined");	
//    }
//    // allocate space for image data and...
//    image_data = (float*) calloc(4 * height * width, sizeof(float));
//    
//    // pull down image from ONLY 1 context
//    display->setCurrent(0);
//    // bind texture
//    glBindTexture(GL_TEXTURE_2D, *(texture_maps->getElement(display->getCurrentContextIndex())));
//    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_FLOAT, image_data);
//    // unbind texture
//    glBindTexture(GL_TEXTURE_2D, 0);
//    
//    for (int x = 0; x < width; x++) {
//        for (int y = 0; y < height; y++) {
//            std::cout << image_data[x+y*width];
//        }
//    }
//    std::cout << "\n";
    
    
	// TODO: update to work with lists
	// TODO: this is wrong, because texture_map is unsigned...
	/*if(texture_map > 0) {
     loaded = true;
     } else {
     loaded = false;
     }*/
}

//void mwMaskStimulus::drawInUnitSquare(StimulusDisplay *display) {
//    //super.drawInUnitSquare(display);
//    std::cout << "drawing in unit square for mask\n";
//    ImageStimulus::drawInUnitSquare(display);
//    std::cout << "finished drawin in unit square for mask\n";
//}

Data mwMaskStimulus::getCurrentAnnounceDrawData() {
    
    //mprintf("getting announce DRAW data for image stimulus %s",tag );
    
    Data announceData(M_DICTIONARY, 9);
    announceData.addElement(STIM_NAME,tag);        // char
    announceData.addElement(STIM_ACTION,STIM_ACTION_DRAW);
    announceData.addElement(STIM_TYPE,STIM_TYPE_IMAGE);
    announceData.addElement(STIM_FILENAME,filename);  
    announceData.addElement(STIM_POSX,last_posx);  
    announceData.addElement(STIM_POSY,last_posy);  
    announceData.addElement(STIM_SIZEX,last_sizex);  
    announceData.addElement(STIM_SIZEY,last_sizey);  
    announceData.addElement(STIM_ROT,last_rot);
    //TODO    announceData.addElement(STIM_ALPHA,last_alpha);  
    
    return (announceData);
}