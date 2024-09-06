addpath('Auxiliary_Toolboxes/KTHorb/ORBIT/Classes')
orbAnal = OrbitDataAnalyser();

filterOutputFolderPath = 'ws/champ_2008_319/EKF/outputs';
cellFilterOutputFile={'20031029.txt','20031029_R.txt'};
cellFilterOutputMatFile={'20031029.mat','20031029_R.mat'};

precise_orbit_folder_path='ws/champ_2003_297_307/UKF/inputs/pod';
cellPrecOrbitDataFile={'2003-10-29.champjpl'};

ref_epoch=[2003 10 29 0 0 10];

measurement_plot_fromDT = [ 2003 10 29 0 0 0 ];
measurement_plot_endDT  = [ 2003 10 29 23 59 50 ];

%-------------------------------------------------------------------
% COMPARE WITH PRECISSE ORBIT
%-------------------------------------------------------------------
orbAnal.readFilterOutputFiles(filterOutputFolderPath,cellFilterOutputFile);
 [absMeanRTN,absRmsRTN,meanRTN,rmsRTN,absMeanXYZ,absRmsXYZ,meanXYZ,rmsXYZ]= ...
 orbAnal.compareToJplPreciseOrbit( precise_orbit_folder_path,cellPrecOrbitDataFile,...
                                   1,[0 0 0],ref_epoch,0)
                               
%-------------------------------------------------------------------
% PLOT ESTIMATED STATE VECTOR PARAMETERS
%-------------------------------------------------------------------
orbAnal.plotstateVectorParam(1,[2003 10 29 0 0 0])

%-------------------------------------------------------------------
% PLOT OBSERVATIONS
%-------------------------------------------------------------------
orbAnal.readFilterMatOutputFile(filterOutputFolderPath,cellFilterOutputMatFile);
orbAnal.plotObservations(measurement_plot_fromDT, measurement_plot_endDT,ref_epoch)
