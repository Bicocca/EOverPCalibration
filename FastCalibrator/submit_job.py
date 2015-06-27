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

name = ["EB"];
#split = ["split","noSplit"];
split = ["noSplit","split"];
cut = ["005"];

for j in range(len(split)):
  for i in range(len(name)):
    print split[j]
    print name[i]
    fn = "Job/Job_"+name[i]+cut[i]+"_"+split[j];
    outScript = open(fn+".sh","w");
    command = "./bin/FastCalibrator"+name[i]+".exe cfg/FastCalibrator_"+name[i]+"_"+split[j]+".cfg";
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'cd '+CMSSWDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+"unbuffer "+command+" > "+fn+"_Norm005.txt");
    outScript.close();
    os.system("chmod 777 "+currentDir+"/"+fn+".sh");
    os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
                                                
