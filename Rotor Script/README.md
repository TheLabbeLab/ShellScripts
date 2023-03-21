**Rotor Scan Input Python Script Tutorial** 

***Motivation***: Typically, to run a Gaussian scan job, you must identify a dihedral of interest on a z-matrix-format input file (.gjf) and write *s 36 10.0* next to it (in addition to changing the route section with the appropriate directions). This can be time-consuming if working with multiple files and increases the likelihood of making a silly mistake somewhere. The goal of these scripts is to identify and generate input (.gjf) files for every rotor of a given molecule (on z- matrix format). It works as well for multiple molecules at the same time. 

There are 3 scripts that you will be using for this: 

- get\_rotors.py 
- zmat2xyz.py 
- get\_geometry.py 

***Step 1:*** Place all 3 .py scripts in the same directory as the molecule you want to generate input files for. For this demo, we will be using Dimethoxymethane.  

D8 ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.001.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.002.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.003.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.004.png)

D1 ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.005.png)

D3 ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.006.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.004.png)

D2 ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.005.png)

Dimethoxymethane: Has 4 rotors. For the numbering above, they will correspond to D3, D1, D2 and D8.   

![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.007.jpeg)

***Step 2:*** Navigate to this directory in your terminal (make sure you have python 3 installed), and run the scripts by entering the following command *python get\_rotors.py* 

![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.008.jpeg)

***Step 3:*** Specify if the file/molecule in question is a Transition State or not (y/n) 

![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.009.jpeg)

***Step 4:*** Specify how many procesors you want to allocate for all the rotor scans for that specific file/molecule. 

![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.010.jpeg)

***Step 5:*** Output! 

Molecule’s  Molecule’s bond  ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.011.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.012.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.013.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.014.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.015.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.016.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.017.png)rotors  

Molecule’ z-matrix  ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.018.png)![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.019.png)

Input files are inside newly created directory. ![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.020.png)

![](Aspose.Words.dd2d1beb-19dc-46c8-9801-f00ca2830e02.021.png)

This process can be done with multiple files at a time, in which case multiple directories will be created with each molecules respective rotor input files inside each one. 
