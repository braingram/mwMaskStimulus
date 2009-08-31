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
        for (int i = 0; i < 4; i++) {
            free(channel_modulus[i]);
        }
        // !!! TODO !!! free fftwf data
    }
}

shared_ptr<Variable> mwMaskStimulus::getRandomSeed() {
    return random_seed;
}

void mwMaskStimulus::makeMask(StimulusDisplay *display) {
    if(!loaded){
        mprintf("Attempted to make mask when image was not loaded, no mask made");
        return;
    }
//    mprintf("-----------------------------------------------------------------");
//    mprintf("---------------------processing %s -------------",tag.c_str());
    
    
//    for (int i = 0; i < 3; i++) {
//        float imageMax = -1000.0;
//        float imageMin = 1000.0;
//        float imageMean = 0.0;
//        for (int j = 0; j < (height * width); j++) {
//            imageMean += image_data[(j*4)+i];
//            imageMax = (image_data[(j*4)+i] > imageMax) ? image_data[(j*4)+i] : imageMax;
//            imageMin = (image_data[(j*4)+i] < imageMin) ? image_data[(j*4)+i] : imageMin;
//        }
//        imageMean /= height * width;
//        mprintf("%i::Pre Mask::::",i);
//        mprintf("\tMax:  %f",imageMax);
//        mprintf("\tMin:  %f",imageMin);
//        mprintf("\tMean: %f",imageMean);
//    }
    
    for (int i = 0; i < 3; i++) {
        float max = -1000.0;
        float min = 1000.0;
        float mean = 0.0;
        for (int j = 0; j < ((height * width)/2+1); j++) {
            mean += channel_modulus[i][j];
            max = (channel_modulus[i][j] > max)? channel_modulus[i][j]: max;
            min = (channel_modulus[i][j] < min)? channel_modulus[i][j]: min;
        }
        mean /= height * width;
        mprintf("%i::Pre Mask Channel Modulus:: %s ::",i,tag.c_str());
        mprintf("%i: Max: %f\tMin: %f\tMean: %f",i,max,min,mean);
    }
    
    
    // ======================================
    //      Making mask fubars colorLines
    // ======================================
    
    
    // generate random phase
    boost::mt19937 rng;
    // !!! TODO !!! should I use an existing pi value?
    boost::uniform_real<> phase_distribution(-3.14,3.14);
    //boost::uniform_real<> phase_distribution(0.0,6.2831853071795862);
    // make generator, which when called, generates a random number between 0 and 1.0
    boost::variate_generator<boost::mt19937, boost::uniform_real<> > random_phase_gen(rng,phase_distribution);
    // allocate space for random phase !!! possibly move this to load !!!
    fftwf_complex *random_phase = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * height * width);
    // use the generator to fill the phase array
    for (int i = 0; i < (height * width); i++) {
        complex<float> temp_phase(0.0,random_phase_gen());
        temp_phase = exp(temp_phase);
        random_phase[i][0] = real(temp_phase);
        random_phase[i][1] = imag(temp_phase);
        //random_phase[i] = reinterpret_cast<fftwf_complex>(exp(temp_phase));
    }
    // make space for mask
    float *mask_data = (float*)calloc(height*width*4,sizeof(float));
    // combine with modulus and do inverse
    for (int i = 0; i < 3; i++) {
        
        // make space for fftw input !!! possibly move this out of the loop? !!!
        //fftwf_complex *in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * height * width);
        fftwf_complex *in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((height * width)/ 2 + 1));
        // fill input
        //for (int j = 0; j < (height * width); j++) {
        for (int j = 0; j < ((height * width)/2 + 1); j++) {
            complex<float> temp_modulus(channel_modulus[i][j],0.0);
            complex<float> temp_phase(random_phase[j][0],random_phase[j][1]);
            complex<float> temp_result;
            temp_result = temp_modulus * temp_phase;
            //in[j] = random_phase[j] * temp_modulus;
            in[j][0] = real(temp_result);
            in[j][1] = imag(temp_result);
        }
        // make space for fftw output !!! possibly move this out of the loop? !!!
        float *out = (float*) calloc(height*width,sizeof(float));
        // make fftw plan !!! possibly move this out of the loop? !!!
        fftwf_plan fft_mask_plan = fftwf_plan_dft_c2r_2d(height, width, in, out, FFTW_MEASURE);
        // exectue fft
        fftwf_execute(fft_mask_plan);
        // copy out to mask_data
        for (int j = 0; j < (height * width); j++) {
            //mask_data[(j*4)+i] = out[j];
            mask_data[(j*4)+i] = out[j]/(height * width);
        }
        
        float max = -1000.0; float min = 1000.0; float mean = 0.0;
        for (int j = 0; j < (height * width); j++) {
            mean += out[j];
            max = (out[j] > max)? out[j]: max;
            min = (out[j] < min)? out[j]: min;
        }
        mean /= (height * width);
        mprintf(" %i :: ifft output :: %s",i,tag.c_str());
        mprintf("    max: %f\tmin: %f\tmean: %f",max,min,mean);
        
        fftwf_destroy_plan(fft_mask_plan);
        fftwf_free(in);
        free(out);
        //lock->unlock();
    }
    
    for (int i = 0; i < 3; i++) {
        float maskMax = -1000.0;
        float maskMin = 1000.0;
        float maskMean = 0.0;
        for (int j = 0; j < (height * width); j++) {
            maskMean += mask_data[(j*4)+i];
            maskMax = (mask_data[(j*4)+i] > maskMax) ? mask_data[(j*4)+i] : maskMax;
            maskMin = (mask_data[(j*4)+i] < maskMin) ? mask_data[(j*4)+i] : maskMin;
        }
        maskMean /= height * width;
        mprintf("%i:: Mask :: %s ::",i,tag.c_str());
        mprintf("%i: Max: %f\tMin: %f\tMean: %f",i,maskMax,maskMin,maskMean);
    }
    
    
    // normalize
    for (int i = 0; i < 3; i++) {
        float maskMax = -1000.0;
        float maskMin = 1000.0;
        float maskMean = 0.0;
        float imageMax = -1000.0;
        float imageMin = 1000.0;
        float imageMean = 0.0;
        for (int j = 0; j < (height * width); j++) {
            maskMean += mask_data[(j*4)+i];
            maskMax = (mask_data[(j*4)+i] > maskMax) ? mask_data[(j*4)+i] : maskMax;
            maskMin = (mask_data[(j*4)+i] < maskMin) ? mask_data[(j*4)+i] : maskMin;
            imageMean += image_data[(j*4)+i];
            imageMax = (image_data[(j*4)+i] > imageMax) ? image_data[(j*4)+i] : imageMax;
            imageMin = (image_data[(j*4)+i] < imageMin) ? image_data[(j*4)+i] : imageMin;
        }
        maskMean /= height * width;
        imageMean /= height * width;
//        mprintf("%i::Mask::::", i);
//        mprintf("\tMax:  %f",maskMax);
//        mprintf("\tMin:  %f",maskMin);
//        mprintf("\tMean: %f",maskMean);
//        mprintf("%i::Image::::",i);
//        mprintf("\tMax:  %f",imageMax);
//        mprintf("\tMin:  %f",imageMin);
//        mprintf("\tMean: %f",imageMean);
        // normalize
        float scale = imageMean/(maskMean-maskMin);
        scale = (scale > ((1-imageMean)/(maskMax-maskMean))) ? ((1-imageMean)/(maskMax-maskMean)) : scale;
        //mScale = min(mean(image)/(mean(mask)-mask.min()), (1-mean(image))/(mask.max()-mean(mask)))
        for (int j = 0; j < (height * width); j++) {
            mask_data[(j*4)+i] = scale * mask_data[(j*4)+i] + (imageMean - maskMean * scale);
        }
        //return mScale*mask+(mean(image)-mean(mask)*mScale)
    }
    
    // copy over alpha
    for (int i = 0; i < (height * width); i++) {
        // -- show original ---
//                mask_data[(i*4)+0] = image_data[(i*4)+0];
//                mask_data[(i*4)+1] = image_data[(i*4)+1];
//                mask_data[(i*4)+2] = image_data[(i*4)+2];
//                mask_data[(i*4)+3] = image_data[(i*4)+3];
        // -- red square --
//        mask_data[(i*4)+0] = 1.0;
//        mask_data[(i*4)+1] = 0.0;
//        mask_data[(i*4)+2] = 0.0;
//        mask_data[(i*4)+3] = 0.5; // 0.0 should be opaque
//        // -- modulus --
//        mask_data[(i*4)+0] = channel_modulus[0][i];
//        mask_data[(i*4)+1] = channel_modulus[1][i];
//        mask_data[(i*4)+2] = channel_modulus[2][i];
        // -- copy alpha --
        mask_data[(i*4)+3] = channel_modulus[3][i];
    }
    
//    // --------- DEBUGGING -------------
//    // print out image statistics
//    float maskMax = 0.0;
//    float maskMin = 1.0;
//    float imageMax = 0.0;
//    float imageMin = 1.0;
//    for (int i = 0; i < (height * width * 4); i++) {
//        maskMax = (mask_data[i] > maskMax)? mask_data[i]: maskMax;
//        maskMin = (mask_data[i] < maskMin)? mask_data[i]: maskMin;
//        imageMax = (image_data[i] > imageMax)? image_data[i]: imageMax;
//        imageMin = (image_data[i] < imageMin)? image_data[i]: imageMin;
//    }
//    mprintf("::::Mask:::: Post Norm");
//    mprintf("\tMax: %f",maskMax);
//    mprintf("\tMin: %f",maskMin);
//    mprintf("::::Image:::: Post Norm");
//    mprintf("\tMax: %f",imageMax);
//    mprintf("\tMin: %f",imageMin);
//    // ---------------------------------
    
    // delete old textures
    // move 'masks' (original images right now) to gpu
    for(int i = 0; i < display->getNContexts(); i++){
		display->setCurrent(i);
        GLuint texture_map;
        texture_map = *(texture_maps->getElement(i));
        //texture_maps.getElement(i,&texture_map);
        // delete old texture
        glDeleteTextures(1,&texture_map);
        // replace with new texture
        glGenTextures(1,&texture_map);
        glBindTexture(GL_TEXTURE_2D, texture_map);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, mask_data);
        gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGBA,width,height,GL_RGBA,GL_FLOAT,mask_data);
		if(texture_map){
			mprintf("Image Mask loaded into texture_map %d", texture_map);
		}
        // !!! do I need to readd the texture to the list !!!
        texture_maps->addElement(i, texture_map);
        glBindTexture(GL_TEXTURE_2D, 0);
	}
    
    // replace with new textures
    fftwf_free(random_phase);
    free(mask_data);
    
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
//    // copy new mask to textures
//    for(int i = 0; i < display->getNContexts(); i++){
//        display->setCurrent(i);
//        texture_maps->addElement(i, texture_map);
//        
//    };
}

void mwMaskStimulus::load(StimulusDisplay *display) {
    // loads image into a texture of format RGBA FLOAT
    
    //GlobalOpenGLContextManager->setGlobalContextCurrent();
    if(loaded){
        // !!! TODO !!!
        // regenerate mask
        makeMask(display);
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
    // BJG: the image data gets loaded upside down (this is probably mac specific)
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
    ilEnable(IL_ORIGIN_SET);
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
    
    // break up image data into R G B A
    for (int i = 0; i < 4; i++) {
        channel_modulus[i] = (float*) calloc(height * width, sizeof(float));
    }
    for (int i = 0; i < (height * width); i++) {
        channel_modulus[0][i] = image_data[(i*4)]; // R
        channel_modulus[1][i] = image_data[(i*4)+1]; // G
        channel_modulus[2][i] = image_data[(i*4)+2]; // B
        channel_modulus[3][i] = image_data[(i*4)+3];// A
    }
    
    for (int i = 0; i < 3; i++) {
        float max = -1000.0;
        float min = 1000.0;
        float mean = 0.0;
        for (int j = 0; j < (height * width); j++) {
            mean += channel_modulus[i][j];
            max = (channel_modulus[i][j] > max)? channel_modulus[i][j]: max;
            min = (channel_modulus[i][j] < min)? channel_modulus[i][j]: min;
        }
        mean /= height * width;
        mprintf("%i:: Pre Modulus :: %s ::",i,tag.c_str());
        mprintf("%i: Max: %f\tMin: %f\tMean: %f",i,max,min,mean);
    }
    
    // ======================================
    //      Modulus fubars thumbs_up
    // ======================================
    
    // !!! TODO !!!
    // Data are in format RGBA and of type FLOAT. So, to put it in a form fftw
    // can understand we must...
    
    // fourier transform image data
    //  plan fft from float -> complex
    //fftwf_plan fft_mask_plan = fftwf_plan_dft_r2c_2d()
    //  plan ifft from complex -> float
    //ifft_mask_plan = fftwf_plan_dft_c2r_2d()
    // store magnitude (not phase)
    for (int i = 0; i < 3; i++) {
        // make space for fftw result
        fftwf_complex *channel_fft = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * height * width);
        //channel_fft[i] = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * height * width);
        // plan fft for this channel
        fftwf_plan fft_mask_plan = fftwf_plan_dft_r2c_2d(height, width, channel_modulus[i], channel_fft, FFTW_MEASURE);
        // need to do some padding
        
        // perform fft storing output in channel_fft[i]
        fftwf_execute(fft_mask_plan);
        // calculate modulus of fft: sqrt(real ** 2 + imag ** 2)
        for (int j = 0; j < (height * width); j++) {
            channel_modulus[i][j] = sqrt(channel_fft[j][0] * channel_fft[j][0] + channel_fft[j][1] * channel_fft[j][1]);
        }
        // !!! TODO !!! do I need scaling? 1/sqrt(rows * cols)
        // clean up
        fftwf_destroy_plan(fft_mask_plan);
        fftwf_free(channel_fft);
        //lock->unlock();
    }
    
    for (int i = 0; i < 3; i++) {
        float max = -1000.0;
        float min = 1000.0;
        float mean = 0.0;
        for (int j = 0; j < ((height * width)/2+1); j++) {
            mean += channel_modulus[i][j];
            max = (channel_modulus[i][j] > max)? channel_modulus[i][j]: max;
            min = (channel_modulus[i][j] < min)? channel_modulus[i][j]: min;
        }
        mean /= height * width;
        mprintf("%i::Post Modulus:: %s ::",i,tag.c_str());
        mprintf("%i: Max: %f\tMin: %f\tMean: %f",i,max,min,mean);
    }
    
    // move 'masks' (original images right now) to gpu
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    for(int i = 0; i < display->getNContexts(); i++){
		display->setCurrent(i);
        
        glGenTextures(1,&texture_map);
        glBindTexture(GL_TEXTURE_2D, texture_map);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, mask_data);
        gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGBA,width,height,GL_RGBA,GL_FLOAT,image_data);
		if(texture_map){
			mprintf("Image Mask loaded into texture_map %d", texture_map);
		}
        // !!! do I need to readd the texture to the list !!!
        texture_maps->addElement(i, texture_map);
        glBindTexture(GL_TEXTURE_2D, 0);
        
//		//texture_map = ilutGLBindMipmaps();
//        texture_maps->addElement(i, texture_map);
//		if(texture_map){
//			mprintf("Mask Image loaded into texture_map %d", texture_map);
//		}
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
    makeMask(display);
    
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