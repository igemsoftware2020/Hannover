#coding latin-1
#********************************************************************************************
#THR 
#semiautodeveloped code
#********************************************************************************************


#<<<<<<< redacted >>>>>>>


#********************************************************************************************
# imports

import cv2
import numpy as np
import random
import os
import math

# custom libraries
from src.bacteria import bacterium


#********************************************************************************************
#general dependency functions
def norm(value):
    if value<=1.0:
        if value>=0.0:
            return value
        else:
            return -value
    else:
        return 1.0


    
    
#********************************************************************************************
#core function
frameDim = (500,750)




def output_directory():
    #get path to python code to retrieve parent folder
    #save to output directory
    path = os.path.dirname(os.path.realpath(__file__))
    parentDir = "\\".join(path.split("\\")[:-1])
    outputDir = parentDir+"\\output\\"
    print("outputDir" + outputDir)
    return outputDir


def storage_file_init(path_out):
    #use corresponding video file path
    #by changing .avi to .txt
    #path_out = video_path.replace("avi", "txt")
    #create simple text file
    txt= open(path_out,"r+")    
    return txt

def storage_file_read_frames(txt):
    return txt.read().split("frame ")
    

def storage_get_bacterium(bacterium_text_line):
    #Check, if the text is likely to contain data
    if("," in bacterium_text_line) and (";" in bacterium_text_line):
        
        bacterium_text_line = bacterium_text_line.replace(" ", "")
        property = bacterium_text_line.split(";")
        #print(property)
        
        #Float conversion
        #index                 = int(float(property[0]))
        totalForce_equivalent = float(property[1])
        width                 = float(property[4])
        length                = float(property[5])
        
        #Boolean conversion
        living = ("true" in property[2].lower())
        moving = ("true" in property[3].lower())

        #Array conversion
        pos              = np.array(property[6].replace("[", "").replace("]", "").split(",")).astype(np.float)
        velocity         = np.array(property[7].replace("[", "").replace("]", "").split(",")).astype(np.float)
        angle            = np.array(property[8].replace("[", "").replace("]", "").split(",")).astype(np.float)
        #velocity_angular = np.array(property[9].replace("[", "").replace("]", "").split(",")).astype(np.float)
        
        return bacterium(pos, width, length, velocity, angle, totalForce_equivalent, living, moving)

#********************************************************************************************
def coreFunction():
    
    #Initialization
    dir = output_directory()
    txt = storage_file_init(dir+"core_test1_1_4.txt")
    
    #Load stored biofilm_data
    frame_data = storage_file_read_frames(txt)
    #store result as new video
    #(out, path_out) = video_init()
    coreLoop(frame_data)

    
        
        
def coreLoop(frame_data):#biofilm, out, txt):  
    #framenumber used for the data storage file
    num_frame = 1
    
    #Loop until the frame is closed
    cv2.imshow("frame", 5) 
    #This unbeatiful call of the frame has to be done
    #In order to exit on exit button the way below
    while(cv2.getWindowProperty("frame", 0) >= 0) and (num_frame<len(frame_data)-1):  
          
    #Draw the frame        
    #Dark_simmulation
        frame = np.zeros((frameDim[0],frameDim[1],3))
        frame[:] = (0.0, 0.25, 0.1)
        
        
         
        lines = frame_data[num_frame].split("\n")
        bacteria = []
        if(len(lines)>2):
            for line_index in range(len(lines)-3):        
                line = lines[1+line_index]
                #print(line)
                bacteria.append(storage_get_bacterium(line))
                
        #Sort bacteria by depth  
        #Redundant later in OpenGL      
        bacteria = sorted(bacteria, key=lambda x: x.pos[2], reverse=True)
        
        bacterium_index = 0
        for Bacterium in bacteria:
            bacterium_index = bacterium_index+1
            #Color coding of the overall-force-sum
            positions = Bacterium.getPositions()
            for pos_index in range(len(positions)):               
                if(pos_index<len(positions)):
                    pos = positions[pos_index]
                else:
                    pos = positions[0] #Draw the zero_circle at last for correct bounding lines
                #Bacteria_height coded by cell-brightness
                z_coding = norm(100.0/(100.0+50-pos[2]+Bacterium.width+0.001))
                #*************************************************************************
                #Front-View
                pos_0 = 0+pos[0]
                pos_1 = 0+pos[1]
                _angle = Bacterium.angle[0]*180/3.14
                
                rotated = True
                if(math.cos(Bacterium.angle[1])<0):
                    rotated = False
              
                if(Bacterium.living==True):
                    b = norm(255*(1.0-1.0*750.0/(750.0+Bacterium.totalForce_equivalent)))
                    frame = cv2.circle(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), int((Bacterium.width*2)), (1.0,norm((0.4+b*0.4+0.1*z_coding)),norm(0.2+0.2*z_coding)), -1)          
                else:
                    b = 1*(1.0-1.0*250.0/(250.0+Bacterium.totalForce_equivalent))
                    frame = cv2.circle(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), int((Bacterium.width*2)), (0.3,(0.25+b),(0.25+b)), -1) 
                
                contour_color = (0.5,0.4,0.2)
                frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                        0, -20-_angle, 20-_angle, contour_color, 2)
                frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                        0, 180-20-_angle, 180+20-_angle, contour_color, 2)
                if((pos_index == len(positions)-1)):
                    if((rotated==False)):
                        frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                                0, 20-_angle, -180+20-_angle, contour_color, 2)
                    else:
                        frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                                0, 180-20-_angle, 20-_angle, contour_color, 2)
                                        
                if((pos_index == 0)):
                    if((rotated==False)):
                        frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                                0, 180-20-_angle, 20-_angle, contour_color, 2)
                    else: 
                        frame = cv2.ellipse(frame, (int(frameDim[1]*0.5+pos_0), int(frameDim[0]*0.5+pos_1)), (int((Bacterium.width*2)),int((Bacterium.width*2))), 
                                                0, 20-_angle, -180+20-_angle, contour_color, 2)    
                             
            #storage_file_append_bacterium(txt, Bacterium, bacterium_index)
                                  
                
        cv2.imshow("frame", frame)
        
        frame = cv2.normalize(frame, None, alpha = 0, beta = 255, norm_type = cv2.NORM_MINMAX, dtype = cv2.CV_32F)
        frame = frame.astype(np.uint8)
        
        #Store frame and data in corresponding files
        num_frame = num_frame + 1
        #out.write(frame)
        
        key = cv2.waitKey(2)
        if(key!=-1):
            if(key==32):
                continue
            print(key)
            break
            #out.release()
            #txt.close()




#********************************************************************************************
#main-method to start the program
#********************************************************************************************
if __name__=="__main__":
    print("started")
    coreFunction()
