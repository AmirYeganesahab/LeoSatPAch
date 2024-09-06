clc
fclose('all');
addpath('Auxiliary_Toolboxes/KTHorb/ORBIT/Classes')
%-----------------------------------------------------------------
% FILTER SETTINGS
%------------------------------------------------------------------
  % Filter workspace for input output
    %filtSet.ws = 'ws/champ_2004_001_031/UKF/';
      % Filter execution time
    filtSet.ws = 'ws/champ_2003_297_307/UKF/';
    filtSet.beginFateTime = FateTime(2003,10,29,0,0,0);
    filtSet.endFateTime   = FateTime(2003,10,30,0,0,0);
 
  % Select the Filter Type
    filtSet.filterType='UnscentedKalman';   
        % Set filterType to 'Kalman' or 'UnscentedKalman'    
    
  % Select the observation type
    filtSet.obsType='Code';   
        % Set filterType to 'Graphic' or 'Code' or 'NavSol'  
    
  % Set State Vector paremeters that will be estimated
    % To estimate atmospheric drag, solar radiation presure, empirical
    % accelertaions and markov process corelletion time
      filtSet.enableDynModParam=1;  
        % Set enableDynModParam to "0" or "1" or "2" 
        % 0 : do not estimate
        % 1 : estimate atmospheric drag coefficient, solar radiation presure
        %     coefficient and empirical accelertaions
        % 2 : estimate atmospheric drag coefficient, solar radiation presure
        %     coefficient, empirical accelertaions and markov process corelletion 
        %     time
        
    % Set initial dynamical model parameters. If estimation of dynamical
    % model parameters are canceled the given values of atmospheric drag  
    % coefficient and solar radiation coefficient will be used as constant 
    % in the filter
  filtSet.initVal.atmDragCoef=2;
  filtSet.initVal.solarRadCoef=1.5;
  filtSet.initVal.empAccel=[1000e-9; 1000e-9; 1000e-9]; % in RTN
  filtSet.initVal.corelTime=600; % in second 
  % Set Satellite Parameters
   filtSet.Satellite.a2m=(1.22/522);  % m2/kg
    filtSet.Satellite.mass=522;       % kg
    filtSet.Satellite.incId=4;
    filtSet.Satellite.Gravity.nmax=70;  
    
  % Set Time System Parameters
    filtSet.TimeRefSystem.UT1_UTC=-0.3652860;     % in second                  
    filtSet.TimeRefSystem.TAI_UTC=32;             % in second
    filtSet.TimeRefSystem.PolarMotion.xp=0.220270;% in second 
    filtSet.TimeRefSystem.PolarMotion.yp=0.242220;% in second   
  
  % Set filter time propagation parameters
    filtSet.Propagation.stepSize=30; % Filter output interval in second
    filtSet.Propagation.timeUpdateThreshold=1; % Update threshold for the  
                                               % propagation step size in second
  % Set statistical parameters of measurements and auxiliary parameters
    filtSet.stat.std.measurementNoise=0.1;   
    filtSet.stat.std.obsSISRE=1.5;  
    
  % Set Initial statistical parameters
    filtSet.stat.std.init.posXYZ=[5 5 5]'; 
    filtSet.stat.std.init.velXYZ=[.05 .05 .05]';      
    filtSet.stat.std.init.atmDragCoeff=0.001;  
    filtSet.stat.std.init.solarRadCoeff=0.001; 
    filtSet.stat.std.init.empAccellRTN=[1e-9  1e-9  1e-9 ]';%[1000e-9 1000e-9 1000e-9]; % m/s2
    filtSet.stat.std.init.corelTime=10;     
    filtSet.stat.std.init.recClcBias=100; 
    filtSet.stat.std.init.ambiguityBias=30;
    
  % Set initialization statistical parameters    
    filtSet.stat.std.initialization.ambiguityBias=25;%50;

  % Set statistical parameters of process noise
    filtSet.stat.std.procNoise.posXYZ=[0 0 0];%[.5 .5 .5]'; % sigma per interval
    filtSet.stat.std.procNoise.velXYZ=[0 0 0];%[.0005 .0005 .005]';       % sigma per interval 
    filtSet.stat.std.procNoise.atmDragCoeff=0.0001;     % sigma for derivative
    filtSet.stat.std.procNoise.solarRadCoeff=0.0001;    % sigma for derivative
    filtSet.stat.std.procNoise.empAccellRTN=[15e-9  15e-9  25e-9 ]; % [30e-9  30e-9  50e-9 ]';%[1000e-9 1000e-9 1000e-9]; % sigma for derivative
    filtSet.stat.std.procNoise.corelTime=0.01;    % sigma for derivative
    filtSet.stat.std.procNoise.recClcBias=1;%2.5;    % sigma for derivative
    filtSet.stat.std.procNoise.ambiguityBias=.01; % m/sqrt(dt)  sigma for derivative  
    
  % Set parameters required for data editing
    % Set the data editing mode.The value of 1 for filtSet.dataEditing.mode
    % indicates recursive outlier detection. The value of 2 is for robust filtering
      filtSet.dataEditing.mode=1;    
    % if recursive outlier detection mode is selected consider the following settings 
      if filtSet.dataEditing.mode==1
         filtSet.dataEditing.outlierFactor=2;      % Outlier Factor
         filtSet.dataEditing.elevationThreshold=10; % elevation threshold in degree
    % if robust mode is selected consider the following settings          
      elseif  filtSet.dataEditing.mode==2
         % Set the level of significance (los) value for the chi-square
         % distribution used to detect faulty measurements. The los can be 
         % set to one of following values ;
         % los=> (.995 .990 .975 .950 .900 .100 .050 .025 .010 .005 )
           filtSet.dataEditig.chiSquare.los=.95; 
      end
    % Minumum number of observation  
      filtSet.dataEditing.minNumObs=4;          
    % GPS recever Antenna Offset from the base location in radial, alongTrack
    % crossTrack directions 
      filtSet.dataEditing.AntOffFromSatBaseLoc=[-0.4306 -1.488 0];  
  
  % For Unscented Kalman Filter; sigma-vector parameters
    filtSet.UKF.kappa=0; 
    filtSet.UKF.alfa=0.1;
    filtSet.UKF.beta=2;   


  % Filter Output Files
    filtSet.enableSaveFilterOutput=1; 

%-----------------------------------------------------------------
% EXUCUTE SIMULATION
%------------------------------------------------------------------
  % Open GPS broadcast ephemerides file
    orbBrdEphFileObj=OrbBrdEphFile(filtSet);
  % Open GPS observation file  
    orbObsFileObj=OrbObsFile(filtSet,'rinex');                   
    outFileObj=OrbFilterOutputFile(filtSet);
    
  % Run the filter
    filterBaseObj=OrbFilter();
    filterObj=filterBaseObj.getInstance(filtSet,orbObsFileObj, ...
                                        outFileObj,orbBrdEphFileObj);
    filterObj.run()
    
    
%   % Save Figures
    fig = figure(filterObj.innovationPlot)
    saveas(fig,'innovationPlot.fig')
%      
    fig = figure(filterObj.residualPlot)
    saveas(fig,'residualPlot.fig')
%     
    fig = figure(filterObj.svnPlot)
    saveas(fig,'svnPlot.fig')
%     
    fig = figure(filterObj.recClcBiasPlot)
    saveas(fig,'recClcBiasPlot.fig')    
    
   
    
