classdef OrbDynamicsModel < handle
    % version 8
    
    properties 
        satDynObj  % 'SatelliteDynamic' object
        filtSet    % filter settings
        timeAndRefSysObj
    end
    
    methods
        
        function this=OrbDynamicsModel(filtSet)
        % FUNCTION
        %   Contructure function to execute orbital dynamics processes
        % INPUTS
        %   filtSet: filter setings
        %----------------------------------------------------------------   
        
        % Create Satellite Dynamics Object
          this.satDynObj=SatelliteDynamic_v8( ...
               filtSet.Satellite.a2m, ...
               filtSet.Satellite.incId, ...
               filtSet.Satellite.mass, ...
               filtSet.Satellite.Gravity.nmax, ...
               filtSet.TimeRefSystem.TAI_UTC, ...
               filtSet.TimeRefSystem.UT1_UTC, ...
               filtSet.TimeRefSystem.PolarMotion.xp, ...
               filtSet.TimeRefSystem.PolarMotion.yp);
        % Store the filter settings
          this.filtSet=filtSet;
        % Create auxiliary object
          this.timeAndRefSysObj=TimeAndRefSystem();
            
        end
        
        function [xp,Pp]=propagateStateCov(this,initSVObj,dt,dh)
        % FUNCTION
        %   Fucntion propagates both state vector and its covariance matrix
        %   through time
        % INPUTS
        %   initSVObj: instance of 'OrbStateVector' class including state
        %              parameters
        %   dt       : elepsed time from the previous time of epoch
        %   dh       : step size for numerical integration
        % OUTPUTS
        %   xp       : predicted state vector
        %   Pp       : predicted covariance matrix
        %---------------------------------------------------------------- 
        
        % Set the state integration parameters
          TT_GPS=19+32.184; % Difference between terrastrial time and GPS time        
          propType='State&Transition';                   % Type of propagation
          stateModelType=this.filtSet.enableDynModParam; % Type of state vector 
        % Set the state vector
          position=initSVObj.position;
          velocity=initSVObj.velocity;
          if stateModelType==0 
             atmDragCoef=this.filtSet.initVal.atmDragCoef(:);
             solarRadCoef=this.filtSet.initVal.solarRadCoef(:);
             empAccel=this.filtSet.initVal.empAccel(:);
             corelTime=this.filtSet.initVal.corelTime(:);
          elseif stateModelType==1 || stateModelType==2 
             atmDragCoef=initSVObj.atmDragCoef(:);
             solarRadCoef=initSVObj.solarRadCoef(:);
             empAccel=initSVObj.empAccel(:); 
             if stateModelType==2 
                corelTime=initSVObj.corelTime(:);  
             else
              corelTime=this.filtSet.initVal.corelTime(:);
             end
          end
        % Propagate the state and transition matrix to the next epoch
          % Modified Julian date of the initial time (in terrastrial time)
            initDTobj=FateTime(initSVObj.stateFateTime);
            initTTobj=plus(initDTobj,TT_GPS);
            MJD_TT0=get(initTTobj,'MJD'); 
          % Time propagation
            [xp,phi]=this.satDynObj.propagate(propType,...
                          stateModelType,...
                          position,...
                          velocity,...
                          atmDragCoef,...
                          solarRadCoef,...
                          empAccel,...
                          corelTime,...
                          MJD_TT0,...
                          dt,...
                          dh);
        % Extend transition matrix and state vector for measurement model parameters
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic')
             sPhi=length(phi);
             % part for receiver clock bias
               phi(sPhi+1,sPhi+1)=1;
               xp=[xp;initSVObj.recClcBias];
             % Transition part for receiver clock bias ambiguity bias
               if strcmp(this.filtSet.obsType,'Graphic')
                  sB=length(initSVObj.ambBias); % Numeber of ambiguity bias 
                                                % parameters
                  phi=[phi,             zeros(sPhi+1,sB);
                       zeros(sB,sPhi+1),eye(sB)         ]; 
                  xp=[xp;initSVObj.ambBias(:)];  
               end
          end
        % Compute predicted covarainace matrix
          Pp=phi*initSVObj.stateCov*phi';          
            
        end

        function [xp]=propagateState(this,initSVObj,dt,dh)
        % FUNCTION
        %   Fucntion propagates both state vector and its covariance matrix
        %   through time
        % INPUTS
        %   initSVObj: instance of 'OrbStateVector' class including state
        %              parameters
        %   dt       : elepsed time from the previous time of epoch
        %   dh       : step size for numerical integration        
        % OUTPUTS
        %   xp       : predicted state vector
        %   Pp       : predicted covariance matrix
        %---------------------------------------------------------------- 
        
        % Set the state integration parameters
          TT_GPS=19+32.184; % Difference between terrastrial time and GPS time        
          propType='State';                   % Type of propagation
          stateModelType=this.filtSet.enableDynModParam; % Type of state vector 
        % Set the state vector
          position=initSVObj.position;
          velocity=initSVObj.velocity;
          [m,n]=size(position);  % Size of state vector
          if stateModelType==0 
             atmDragCoef=ones(1,n)*this.filtSet.initVal.atmDragCoef;
             solarRadCoef=ones(1,n)*this.filtSet.initVal.solarRadCoef;
             empAccel=repmat(this.filtSet.initVal.empAccel(:),[1,n]);
             corelTime=ones(1,n)*this.filtSet.initVal.corelTime;
          elseif stateModelType==1 || stateModelType==2 
             atmDragCoef=initSVObj.atmDragCoef;
             solarRadCoef=initSVObj.solarRadCoef;
             empAccel=initSVObj.empAccel; 
             if stateModelType==2 
                corelTime=initSVObj.corelTime;  
             else
              corelTime=ones(1,n)*this.filtSet.initVal.corelTime;
             end
          end
        % Propagate the state and transition matrix to the next epoch
          % Modified Julian date of the initial time (in terrastrial time)
            initDTobj=FateTime(initSVObj.stateFateTime);
            initTTobj=plus(initDTobj,TT_GPS);
            MJD_TT0=get(initTTobj,'MJD'); 
          % Time propagation
            for i=1:n
                [xp(:,i)]=this.satDynObj.propagate(propType,...
                          stateModelType,...
                          position(:,i),...
                          velocity(:,i),...
                          atmDragCoef(:,i),...
                          solarRadCoef(:,i),...
                          empAccel(:,i),...
                          corelTime(:,i),...
                          MJD_TT0,...
                          dt,...
                          dh);
            end
        % Extend the state vector for measurement model parameters
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic')
             % part for receiver clock bias
               xp=[xp;initSVObj.recClcBias];
             % Transition part for receiver clock bias ambiguity bias
               if strcmp(this.filtSet.obsType,'Graphic')
                  xp=[xp;initSVObj.ambBias];  
               end
          end         
            
        end        
        
        function [Q]=getProcessNoise(this,initSVObj,dt)
        % FUNCTION
        %   Compute process noise
        % INPUTS
        %   initSVObj: instance of 'OrbStateVector' class including state
        %              parameters
        %   dt       : elepsed time from the previous epoch
        % OUTPUTS
        %   Q        : Computed process noise
        %----------------------------------------------------------------               
        
        % Construct the process noise covariance matrix for position and
        % velocity vector components
          if this.filtSet.enableDynModParam==0
             Qdiag=[ (this.filtSet.stat.std.procNoise.posXYZ(:)').^2, ...
                     (this.filtSet.stat.std.procNoise.velXYZ(:)').^2];
             Q=diag(Qdiag);
          end
        % If filter are set to estimate dynamical model parameters which
        % including empirical acceleration parameters, compute the process
        % noise covariance by considering the dynamic model compensation method
          % Set the parameters
            if this.filtSet.enableDynModParam==1
               corelTime=this.filtSet.initVal.corelTime;
            elseif this.filtSet.enableDynModParam==2
                corelTime=initSVObj.corelTime;
            end
          % Compute the process noise covariance 
            if this.filtSet.enableDynModParam==1 || ...
               this.filtSet.enableDynModParam==2 
               % Process noise covariance for the position, velocity and
               % empirical accelerations
               % Reference: Lee, D. (2005). Nonlinear Bayesian filtering with 
               % applications to estimation and navigation. Texas A&M University. 
                 sigma_emp=(this.filtSet.stat.std.procNoise.empAccellRTN(:)).^2*corelTime/2;
                 L=diag(sigma_emp*(1-exp(-2*dt/corelTime)));
                 Q=zeros(11,11);  % Initiaize the process noise covariance matrix
% %                  Q(1:3,1:3)=dt^4/4*L;    Q(1:3,4:6)=dt^3/2*L; Q(1:3,9:11)=dt^2/2*L;
% %                  Q(4:6,1:3)=Q(1:3,4:6);  Q(4:6,4:6)=dt^2*L;   Q(4:6,9:11)=dt*L;
% %                  Q(9:11,1:3)=Q(1:3,9:11);Q(9:11,4:6)=dt*L;    
                 Q(9:11,9:11)=L;
               % Process Noise covaraince for corellation time
                 if this.filtSet.enableDynModParam==2
                    Q(12,12)=this.filtSet.stat.std.procNoise.corelTime^2*dt;   
                 end
               % Process Noise covaraince for the atmospheric drag and solar
               % radiation pressure coefficients 
                 Q(7,7)=this.filtSet.stat.std.procNoise.atmDragCoeff^2*dt; 
                 Q(8,8)=this.filtSet.stat.std.procNoise.solarRadCoeff^2*dt; 
% %                % Convert the position, velocity components of the covariance  
% %                % matrix from RTN system to the ECEF system
% %                  pos=initSVObj.position; vel=initSVObj.velocity;
% %                  orb2geo=this.timeAndRefSysObj.orb2geo(pos(1),pos(2),pos(3),...
% %                                                  vel(1),vel(2),vel(3)); 
% %                  Q(1:3,1:3)=orb2geo*Q(1:3,1:3)*orb2geo';
% %                  Q(1:3,4:6)=orb2geo*Q(1:3,4:6)*orb2geo'; 
% %                  Q(1:3,9:11)=orb2geo*Q(1:3,9:11)*orb2geo'; 
% %                  
% %                  Q(4:6,1:3)=orb2geo*Q(4:6,1:3)*orb2geo';
% %                  Q(4:6,4:6)=orb2geo*Q(4:6,4:6)*orb2geo'; 
% %                  Q(4:6,9:11)=orb2geo*Q(4:6,9:11)*orb2geo';                  
             
            end
        % Compute process noise covariance for measurement model parameters
          if strcmp(this.filtSet.obsType,'Code') || ...
             strcmp(this.filtSet.obsType,'Graphic') 
             sQ=length(Q); % Size of process noise
             % Process noise component for receiver clock bias
               Q(sQ+1,sQ+1)=this.filtSet.stat.std.procNoise.recClcBias^2*dt;
             % Process noise component for ambiguity bias terms of GRAPHIC 
             % observables  
               if strcmp(this.filtSet.obsType,'Graphic') 
                  [rB,cB]=size(initSVObj.ambBias); % Number of ambiguity bias terms
                  Q(sQ+2:sQ+2+rB-1,sQ+2:sQ+2+rB-1)=eye(rB)* ...
                          this.filtSet.stat.std.procNoise.ambiguityBias^2*dt;
               end
          end

          
      
        end
        
   

        
    end
    
    
end

