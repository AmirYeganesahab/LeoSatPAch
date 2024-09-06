classdef OrbStateVector < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Orbit State Parameters
          position
          velocity
        % Dynamic Model Parameters
          atmDragCoef
          solarRadCoef
          empAccel
        % Adaptive Filter Parameters
          corelTime
        % Measurement Model Parameters
          recClcBias
          ambBias
        % Compound State Vector
          stateVector
        % State Covariance
          stateCov
        % Weight (for only unscneted Kalman filter)
          stateWeight 
        % Covariance scaling matrix
          covSclMat
        % Innovations
          innov
        % Date and time
          stateFateTime
        % Update Flag
          updateFlag 
        % Transformation vector from body system to GPS system at ECEF
        % system
          bodyToGps
    end
    
    methods
        function this=OrbStateVector(varargin)
                                 
        % FUNCTION
        %   Contructure function for the class to store state vector
        %   parameters of the orbit
        % INPUTS 
        %   Function can be called in two way;
        %      - obj=OrbStateVector(position,velocity,atmDragCoef,solarRadCoef, ...
        %                           empAccel,corelTime,recClcBias,ambBias,...
        %                           stateFateTime )
        %      OR
        %      - obj=OrbStateVector(stateVect,enableDynModParam,obsType, ...
        %                           stateFateTime)
        %      WHERE  
        %         position     : position vector
        %         velocity     : velocity vector
        %         atmDragCoef  : atmospheric drag coefficient
        %         solarRadCoef : solar radiation coefficient
        %         empAccel     : empirical acceleration vector
        %         corelTime    : corelation time for Gauss Markov Process
        %         recClcBias   : receiver clock bias
        %         ambBias      : ambiguity bias vectror for frequency measurements 
        %         stateFateTime: array including Date and time, such as;
        %                        [year month day hour minute second]
        %         stateVect    :
        %         obsType      :        
        %         enableDynModParam :
        %----------------------------------------------------------------       
        
        % Set the initial values
          position=[];velocity=[];atmDragCoef=[];solarRadCoef=[];
          empAccel=[];corelTime=[];recClcBias=[];ambBias=[]; 
        % Set the class variables
          n_var=length(varargin);
          if n_var==9
             % Orbit State Parameters
               position=varargin{1};
               velocity=varargin{2};
             % Dynamic Model Parameters
               atmDragCoef=varargin{3};
               solarRadCoef=varargin{4};
               empAccel=varargin{5};
             % Adaptive Filter Parameters
               corelTime=varargin{6};
             % Measurement Model Parameters
               recClcBias=varargin{7};
               ambBias=varargin{8};
             % Date and time
               stateFateTime=varargin{9};     
          elseif n_var==4
              % Get the input variables
                stateVect=varargin{1};
                enableDynModParam=varargin{2};
                obsType=varargin{3};
              % Set position and velocity
                position=stateVect(1:3,:); velocity=stateVect(4:6,:);
                stInd=6;
              % Set dynamical model parameters
                if enableDynModParam==1  || enableDynModParam==2 
                   atmDragCoef=stateVect(7,:); 
                   solarRadCoef=stateVect(8,:);
                   empAccel=stateVect(9:11,:);
                   stInd=11;
                   if enableDynModParam==2
                      corelTime=stateVect(12,:); 
                      stInd=12;
                   end
                end
              % Set measurement model parameters
                if strcmp(obsType,'Code') || strcmp(obsType,'Graphic')
                   recClcBias=stateVect(stInd+1,:); 
                   if strcmp(obsType,'Graphic')
                      ambBias=stateVect(stInd+2:end,:);
                   end
                end
              % Set the date time parameter
                stateFateTime=varargin{4};            
          else
              error('Wrong inputs')
          end
        % Store the class variables
          % Orbit State Parameters
            this.position=position;
            this.velocity=velocity;
          % Dynamic Model Parameters
            this.atmDragCoef=atmDragCoef;
            this.solarRadCoef=solarRadCoef;
            this.empAccel=empAccel;
          % Adaptive Filter Parameters
            this.corelTime=corelTime;
          % Measurement Model Parameters
            this.recClcBias=recClcBias;
            this.ambBias=ambBias;
          % Date and time
            this.stateFateTime=stateFateTime;     
          % Set the compond State Vector
            this.stateVector=[position;velocity;atmDragCoef; ...
                              solarRadCoef;empAccel;corelTime;...
                              recClcBias;ambBias];         
        end
    end
    
end

