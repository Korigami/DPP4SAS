# DPP4SAS
Implementation of algotithms for sampling from Determinantal Point Processes (DPPs) for SAS. 

We provide the implementation in three parts:
* Set of 4GL macros,
* Package in SAS/IML,
* Enterprise Miner nodes.

## Installation of IML package

To install the package use within proc iml:
```
package install 'path to dppsampl.zip';
```
Load the package with:
```
package load dppsampl;
```

## Installation guide for SAS Enterprise Miner

You might need administrator privilages to perform this operation.
1. Clone this repository or download it as ZIP, unpack it, and place the contents in the location on your disk available for SAS. 

2. Open SAS session, run the code from the file DPP4SAS\SASMINER\create_sources3.sas

3. Copy XML files from repository (DPP4SAS\SASMINER\*.xml) to directory SASHome\SASEnterpriseMinerWorkstationConfiguration\\&lt;version&gt;\WEB-INF\classes\components.

Usually SASHome directory is in C:\Program Files directory on Windows.

&lt;version&gt; means wersion of SAS Enterprise Miner, e.g. 15.1

4. Add below line to file SASHome\SASEnterpriseMinerWorkstationConfiguration\\&lt;version&gt;\WEB-INF\classes\components\EMList.txt
```
SpectralClustering=SpectralClustering.xml
```

5. Copy icons from 
 * DPP4SAS\SASMINER\gif, 
 * DPP4SAS\SASMINER\gif16 
 * DPP4SAS\SASMINER\gif32 

to repsective directories in SASHome\SASEnterpriseMinerWorkstationConfiguration\\&lt;version&gt;\WEB-INF\classes\components.

You might need administrator privilages to perform this operation.

6. Start or restart SAS Enterprise Miner

7. In Exploration tab there should be available new Spectral Clustering node

