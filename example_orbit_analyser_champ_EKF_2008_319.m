addpath('Auxiliary_Toolboxes/KTHorb/ORBIT/Classes')
orbAnal = OrbitDataAnalyser();

filterOutputFolderPath = 'ws/champ_2008_319/EKF/outputs';
cellFilterOutputFile={'20081114.txt'};
cellFilterOutputMatFile={'20081114.mat'};

precise_orbit_folder_path='ws/champ_2008_319/EKF/inputs/pod';
cellPrecOrbitDataFile={'2008-11-14.champjpl'};

ref_epoch=[2008 11 14 0 0 10];

measurement_plot_fromDT = [ 2008 11 14 0 0 0 ];
measurement_plot_endDT  = [ 2008 11 14 23 59 50 ];

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
orbAnal.plotstateVectorParam(1,[2008 11 14 0 0 0])

%-------------------------------------------------------------------
% PLOT OBSERVATIONS
%-------------------------------------------------------------------
orbAnal.readFilterMatOutputFile(filterOutputFolderPath,cellFilterOutputMatFile);
orbAnal.plotObservations(measurement_plot_fromDT, measurement_plot_endDT,ref_epoch)
