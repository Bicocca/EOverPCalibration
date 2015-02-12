#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = "/gwpool/users/brianza/PHD/CMSSW_6_1_1/src/";
ReducedTreeDir = "";

name = ["EE"];

for i in range(len(name)):
    fn = "Job/Job_"+name[i];
    outScript = open(fn+".sh","w");
    command = "./bin/FastCalibrator"+name[i]+".exe cfg/FastCalibrator_"+name[i]+"_split.cfg";
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'cd '+CMSSWDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+"unbuffer "+command+" > "+fn+"_output_Oct22_Run2012D_Cal_Dic2012_noOldCorr.txt");
    outScript.close();
    os.system("chmod 777 "+currentDir+"/"+fn+".sh");
    os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
                                                
