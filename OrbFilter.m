classdef OrbFilter < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        filtSet
        obsFileObj
        brdEphFileObj    
        measModelObj
        dynModelObj
        timeRefSysObj
        filterInitTime
        outFileObj
    end
    
    methods
        function this=OrbFilter()
            addpath('Auxiliary_Toolboxes/KTHorb/ORBIT/Classes')
          
        end
        
        function filterObj=getInstance(this,varargin)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Creates object of specified type of filter; 'Kalman Filter', or
        %   'Unscented Kalman Filter'
        % INPUTS
        %   varargin : cell structure including input parameters. The class
        %              can be initiated as follows;
        %
        %              For Code and GRAPHIC measurements
        %                 obj=OrbFilter(filtSet,obsFileObj,outFileObj,brdEphObj)
        %              For navigation solution measurements
        %                 obj=OrbFilter(filtSet,obsFileObj,outFileObj)    
        %              
        %              where
        %                filtSet   : data structure including filter settings
        %                obsFileObj: Handle to instance of 'OrbObsFile' class
        %                brdEphObj: Handle to instance of 'BrdEphObj' class  
        %                outFileObj: Handle to instance of 'OrbFilterOutputFile'
        %                            to store filter outputs
        % OUTPUTS 
        %   ModelObj : Handle to child class instance
        %--------------------------------------------------------------  
        
        % Get and et the filter setings
          filtSet=varargin{1};
        % Create filter Object
          if strcmp(filtSet.filterType,'Kalman')
             filterObj=OrbKalmanFilter();
          elseif strcmp(filtSet.filterType,'UnscentedKalman')
             filterObj=OrbUnscentedKalmanFilter();
          else
             error ('Input filter type are not supported')
          end
        % Store filter settings
          filterObj.filtSet=filtSet;
        % Create Measurement Model Objects
          baseMeasModelObj=OrbMeasurementModel();
          if strcmp(filtSet.obsType,'NavSol')
             filterObj.measModelObj=baseMeasModelObj.getInstance(...
                                           filtSet,varargin{2});
          elseif strcmp(filtSet.obsType,'Graphic') || ...
                 strcmp(filtSet.obsType,'Code')
                 filterObj.measModelObj=baseMeasModelObj.getInstance( ...
                                         filtSet,varargin{2},varargin{4});
          else
             error ('Input measurement model type can not be recognised')
          end
        % Create dynamical model object
          filterObj.dynModelObj=OrbDynamicsModel(filtSet);
        % Set data files
          filterObj.obsFileObj=varargin{2};     
          filterObj.outFileObj=varargin{3};
          n=length(varargin);  
          if n==4
            filterObj.brdEphFileObj=varargin{4};
          end  
        % Set time system object
          filterObj.timeRefSysObj=TimeAndRefSystem();
            
        end
        
        function run(this)
 
        % Set Auxiliary variables and constants    
          dtPropStep=this.filtSet.Propagation.stepSize; % Filter propagation step size
          dtPropThreshold=this.filtSet.Propagation.timeUpdateThreshold; 
          
        % Get initial state vector
          [initSVobj,initEDobj]=this.getInitialStateVector();
        % Date and time of current epoch
          initDTobj=FateTime(initSVobj.stateFateTime);
          this.filterInitTime=initSVobj.stateFateTime;

        % Run the filter loop
          while 1
            % Get observation data at next epoch
              % Skip the observations acquired before the propagation time step
                while 1
                  % Get the observations
                    [nextEDobj,eof_flag]= ...
                          this.obsFileObj.getNextObs(this.filtSet.obsType);
                    % Control end of file
                      if eof_flag; break ; end                        
                  % Date and time of the current observation epoch
                    nextDTobsObj=FateTime(nextEDobj.epochFateTime);                       
                    dtObs=minus(nextDTobsObj,initDTobj);  
                  % Control the observation loop
                    if  dtObs+dtPropThreshold>=dtPropStep; break; end                   
                end
              % Control end of file
                if eof_flag; fclose('all');break ; end  

                

          
            % Propagation of state vector through the time
              do=1;
              while do~=0                 
                % Compute the propagation time
                  if dtObs<=dtPropStep+dtPropThreshold
                      dtPropTime=dtObs;
                      do=0;
                  else
                      dtPropTime=dtPropStep;
                  end   
                % Call the timeUpdate function of the child class for
                % propagation
                  predSVobj=this.timeUpdate(initSVobj,dtPropTime); 
                % If only prediction will be done, set the time and state
                % vector for next loop and save the predicted state
                  if do~=0
                     initSVobj=predSVobj;
                     initDTobj=FateTime(predSVobj.stateFateTime);
                     dtObs=minus(nextDTobsObj,initDTobj);  
                     % Store the predicted filter outputs
                       this.outFileObj.saveStateDataIntoTextFile(initSVobj)                     
                  end

                  
              end
              strcat('Time Update Done  : ',num2str(predSVobj.stateFateTime))              
            % Measurement Update step. Call the 'measurementUpdate' function 
            % of the child class
              [updSVobj,updNextEDobj,updateFlag]=this.measurementUpdate( ...
                                    initEDobj,nextEDobj,predSVobj);  
                                                        %,initSVobj)
                            
                                
            % Set the time and state vector for next loop
              if updateFlag==0
                 initDTobj=FateTime(updSVobj.stateFateTime);                      
                 initSVobj=updSVobj;  
                 initEDobj=updNextEDobj;
                 strcat('Measurement Update Done  : ',num2str(updSVobj.stateFateTime))                 
              else
                 initDTobj=FateTime(predSVobj.stateFateTime);                      
                 initSVobj=predSVobj;  
                 strcat('Only Prediction Done  : ',num2str(predSVobj.stateFateTime))                 
              end
            % Store the filter outputs
              this.outFileObj.saveStateDataIntoTextFile(initSVobj)
              if this.filtSet.enableSaveFilterOutput==1
                 this.outFileObj.saveFilterDataIntoMatFile(initSVobj,updNextEDobj)
              end
            % Check the filter end time
              if get(initDTobj,'MJD')>=get(this.filtSet.endFateTime,'MJD')
                      disp('END OF FILTER')
%                       close('all')
                      break;
              end
              
          end
          
        end
        
        
        function [stateVectorObj,epochDataObj]=getInitialStateVector(this)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Computes the initial state vector.
        % INPUTS
        %    
        % OUTPUTS 
        %   stateVectorObj : instance of 'OrbStateVector' class including 
        %                    initial state vector parameters   
        %   epochDataObj   : instance of 'OrbEpochData' class including
        %                    observations and auxiliary data at initial
        %                    time of epoch
        %--------------------------------------------------------------               
        
        % Set auxiliary variables
          c = 299792458;          % speed of light in vacum (m/sec)
          lambda_L1=(c/1575.42e6); % wavelength of L1 carrier
        % Initialize the state vector parameters
          position=[];velocity=[];
          atmDragCoef=[];solarRadCoef=[]; empAccel=[]; corelTime=[]; 
          recClcBias=[]; ambBias=[]; stateFateTime=[];        
        % Find the initial position and velocity
          % For Code and GRAPHIC observations
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic')
            % Read the data at first three epoch and compute the kinematic
            % position
              i=1;
%               'Filter Initialization' %#ok<NOPRT>
              while i<4
                % Read observation data
                  if strcmp(this.filtSet.obsType,'Code')
                     [orbEDobj{i},health_flag]= ...
                                       this.obsFileObj.getNextObs('Code');
                  elseif strcmp(this.filtSet.obsType,'Graphic')
                     [orbEDobj{i},health_flag]= ...
                                       this.obsFileObj.getNextObs('Graphic');  
                  end
                  if get(FateTime(orbEDobj{i}.epochFateTime),'MJD')<get(this.filtSet.beginFateTime,'MJD')
                      %disp(strcat('Search the observation file for the filter beginning,obs time:', ...
                                   %num2str(orbEDobj{i}.epochFateTime)))
                      continue;
                  end
                % Compute kinematic position
                  if health_flag==0 % Skip the unhealty data
                     [x{i}]=this.kinematicPositioning(orbEDobj{i}.epochFateTime, ...
                                                      orbEDobj{i}.obs.C1,...
                                                      orbEDobj{i}.svn); 
                     i=i+1; 
                  end
              end
              position=x{2}(1:3); stateFateTime=orbEDobj{2}.epochFateTime; 
              epochDataObj=orbEDobj{2};              
          % For Code and GRAPHIC observations              
           elseif strcmp(this.filtSet.obsType,'NavSol') 
                % Read observation data
                  [orbEDobj,health_flag]= this.obsFileObj.getNextObs('NavSol');
                  position=orbEDobj.obs.NavSol(1:3);
                  velocity=orbEDobj.obs.NavSol(4:6);
                  epochDataObj=orbEDobj;    
                  stateFateTime=orbEDobj.epochFateTime; 
           else
               error('Observation type are not supported')              
           end

        % Compute velocity at second epoch by numerical difrentiation
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic')        
             % Convert the time into second and normalized 
               epo1=FateTime(orbEDobj{1}.epochFateTime);
               epo2=FateTime(orbEDobj{2}.epochFateTime);
               epo3=FateTime(orbEDobj{3}.epochFateTime);
               t3s=minus(epo3,epo1);t2s=minus(epo2,epo1);t1s=minus(epo1,epo1);
             % Use numerical difrentiation to compute velocity at middle point
               v=this.numDiffSOLag(t1s,t2s,t3s,x{1}(1:3,1),x{2}(1:3,1),x{3}(1:3,1),t2s); 
               velocity=v; 
          end
        % Set Dynamic Model Parameters
          if (this.filtSet.enableDynModParam==1) || ...
             (this.filtSet.enableDynModParam==2)
             % Set Initial atmospheric drag coefficient
               atmDragCoef=this.filtSet.initVal.atmDragCoef(:);
             % Set Initial solar radiation coefficient
                solarRadCoef=this.filtSet.initVal.solarRadCoef(:);
             % Set empirical accelerations
               empAccel=this.filtSet.initVal.empAccel(:);
             % Set Markov process corelletaion time
              if (this.filtSet.enableDynModParam==2)
                 corelTime=this.filtSet.initVal.corelTime(:);
              end
          end
        % Set Measurement Model Parameters
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic')
             % Set the receiver clock bias
               recClcBias= x{2}(4);
              
             % Compute the initial ambiguity bias for 'Graphic' observables
               if strcmp(this.filtSet.obsType,'Graphic')
                  % (Code range - Phase Range) 
                    ambBias=(orbEDobj{2}.obs.C1-(orbEDobj{2}.obs.L1*lambda_L1))./2;
               end
          end
        % Set the initial covariance matrix
           std_diag=[];
          % Set covariance matrix for position and velocity
            std_diag=[ std_diag; this.filtSet.stat.std.init.posXYZ(:); ...
                                 this.filtSet.stat.std.init.velXYZ(:)];
          % Set covariance matrix for dynamic model parameters
            if (this.filtSet.enableDynModParam==1) || ...
               (this.filtSet.enableDynModParam==2)                    
                std_diag=[ std_diag; this.filtSet.stat.std.init.atmDragCoeff(:); ...
                                     this.filtSet.stat.std.init.solarRadCoeff(:); ...
                                     this.filtSet.stat.std.init.empAccellRTN(:) ];
                if (this.filtSet.enableDynModParam==2)
                   std_diag=[ std_diag; this.filtSet.stat.std.init.corelTime(:)];
                end
            end
          % Set covariance matrix for measurement model parameters
            if strcmp(this.filtSet.obsType,'Code') || ...
               strcmp(this.filtSet.obsType,'Graphic')
               std_diag=[std_diag;this.filtSet.stat.std.init.recClcBias(:)];
               if strcmp(this.filtSet.obsType,'Graphic')
                  std_diag=[std_diag;ones(length(ambBias),1)* ...
                            this.filtSet.stat.std.init.ambiguityBias(:)];
               end
            end
          % Set the whole covariance matrix
            P=diag(std_diag.^2);
        % Create an instance of 'OrbStateVector' class to store variables
          % Store the state parameters   
            stateVectorObj=OrbStateVector(position,velocity,atmDragCoef, ...
                                          solarRadCoef, empAccel,corelTime, ...
                                          recClcBias, ambBias, stateFateTime);                                  
          % Store the covariance matrix
           stateVectorObj.stateCov=P;
        end
        
        function [x]=kinematicPositioning(this,t,z,sv_no)
        %-----------------------------------------------------------------            
        % FUNCTION
        %   Computes position, velocity and receiver clock bias
        % INPUTS
        %   t    : date and time, such as [year month day hour minute second] 
        %   z    : pseudo code observation
        %   sv_no: vehicle number of GPS satellites
        % OUTPUTS 
        %  x     : vector including estimated position and receiver clock 
        %          correction, such ad;    x=[x y z clc_rec]
        %--------------------------------------------------------------                 
            

            x=zeros(4,1);
            for i=1:10
              % Compute satellite position
                [GPS_pos,GPSclc_corr]=this.brdEphFileObj.getGpsSatParam( ...
                                      t,x(1:3),sv_no,x(4));  
                                  
              % Compute geometric distance,d
                zg = sqrt((power(GPS_pos(:,1)-x(1),2)+ ...
                           power(GPS_pos(:,2)-x(2),2)+ ...
                           power(GPS_pos(:,3)-x(3),2)));
              % Measurement partials      
                H = zeros(length(z),length(x));
                for j = 1:length(z)
                    H(j,1) =-(GPS_pos(j,1)-x(1))/zg(j);
                    H(j,2) =-(GPS_pos(j,2)-x(2))/zg(j);
                    H(j,3) =-(GPS_pos(j,3)-x(3))/zg(j);
                    H(j,4) = 1;
                end 
              % Approximate range
                zc=zg+x(4)-GPSclc_corr;% computed range
              % Least Square Estimation
                dz=(z-zc);
                dx=pinv(H'*H)*H'*dz;
                x=x+dx;
            end
        end
        
        function [dx_dt]=numDiffSOLag(this,t1,t2,t3,x1,x2,x3,t)
        % FUNCTION
        %    computes the derivative of tabulated data using difreantiation of
        %    second order lagrange interpolant
        % INPUTS
        %    t1,t2,t3 : time of given data 
        %    x1,x2,x3 : data 
        %    t        : interpolation time
        % OUTPUTS
        %    dx_dt    : time derivative at time t
        %----------------------------------------------------------------
        
          a=x1.*((2*t-t2-t3)/((t1-t2)*(t1-t3)));
          b=x2.*((2*t-t1-t3)/((t2-t1)*(t2-t3)));
          c=x3.*((2*t-t1-t2)/((t3-t1)*(t3-t2)));
          dx_dt=a+b+c;
        end
        
        
 
        
        
    end
    
end

