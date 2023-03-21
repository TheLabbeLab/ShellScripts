#!/bin/bash
#make Tansfer file
mkdir Transfer
mkdir Mess_Inputs
mkdir Log_files

#cd into the Transfer folder and make the Log_files folder 
cp *.mess Mess_Inputs
cp -r me me Mess_Inputs
#copy all the log folders 
cp *.log Log_files

#Transfer the Pesviewer.inp, XYZ folder, and mess files folder into the Transfer folder
cp pesviewer.inp pesviewer.inp Transfer
cp -r xyz xyz Transfer 
cp -r Mess_Inputs Mess_Inputs Transfer
cp -r Log_files Log_files Transfer
rm -r Log_files
rm -r Mess_Inputs



