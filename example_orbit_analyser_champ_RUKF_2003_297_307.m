addpath('Auxiliary_Toolboxes/KTHorb/ORBIT/Classes')
orbAnal = OrbitDataAnalyser();

filterOutputFolderPath = 'ws/champ_2003_297_307/RUKF/outputs';
cellFilterOutputFile={'20031029.txt'};
cellFilterOutputMatFile={'20031029.mat'};

precise_orbit_folder_path='ws/champ_2003_297_307/RUKF/inputs/pod';
cellPrecOrbitDataFile={'2003-10-29.champjpl'};

ref_epoch=[2003 10 29 0 0 10];

covariance_scale_factor_fromDT = [ 2003 10 29 22 35 0  ];
covariance_scale_factor_endDT  = [ 2003 10 29 23 3 20  ];

measurement_plot_fromDT = [ 2003 10 29 0 0 0 ];
measurement_plot_endDT  = [ 2003 10 30 0 0 0 ];

%-------------------------------------------------------------------
% COMPARE WITH PRECISSE ORBIT
%-------------------------------------------------------------------
orbAnal.readFilterOutputFiles(filterOutputFolderPath,cellFilterOutputFile);
 [absMeanRTN,absRmsRTN,meanRTN,rmsRTN,absMeanXYZ,absRmsXYZ,meanXYZ,rmsXYZ]= ...
 orbAnal.compareToJplPreciseOrbit( precise_orbit_folder_path,cellPrecOrbitDataFile,...
                                   1,[0 0 0],ref_epoch,0)
%[absMeanRTN,absRmsRTN,meanRTN,rmsRTN,absMeanXYZ,absRmsXYZ,meanXYZ,rmsXYZ]= ...
% orbAnal.compareToJplPreciseOrbit( precise_orbit_folder_path,cellPrecOrbitDataFile,...
                                %   1,[0 0 0],ref_epoch,0)                               
                               
%-------------------------------------------------------------------
% PLOT ESTIMATED STATE VECTOR PARAMETERS
%-------------------------------------------------------------------
%orbAnal.plotstateVectorParam(1,[2003 10 29 0 0 0])

%-------------------------------------------------------------------
% PLOT MEASUREMENT COVARIANCE SCALE MATRIX
%-------------------------------------------------------------------
orbAnal.readFilterMatOutputFile(filterOutputFolderPath,cellFilterOutputMatFile);
orbAnal.plotCovSclMatrix(covariance_scale_factor_fromDT,covariance_scale_factor_endDT,ref_epoch)

%-------------------------------------------------------------------
% PLOT OBSERVATIONS
%-------------------------------------------------------------------
%orbAnal.readFilterMatOutputFile(filterOutputFolderPath,cellFilterOutputMatFile);
%orbAnal.plotObservations(measurement_plot_fromDT, measurement_plot_endDT,ref_epoch)
