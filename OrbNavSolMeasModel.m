classdef OrbNavSolMeasModel < OrbMeasurementModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function this=OrbNavSolMeasModel()
            
        end
        
        function [updSVobj,updEDobj,minObsFlag]=dataEditing(this,...
                                            currEDobj,predSVobj)
        %-----------------------------------------------------------------                                        
        % FUNCTION
        %   Function search for the observation, and  than removes the outliers
        % INPUTS
        %   currEDobj : instance of 'OrbEpochData' class including observations
        %               GPS ephemerides and other auxiliary data at current
        %               time   
        %   predSVobj : instance of 'OrbStateVector' class including predicted 
        %               state vector parameters  
        % OUTPUTS
        %   updSVobj  : instance of 'OrbStateVector' class including updated 
        %               state vector parameters and covariance matrix           
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %   minObsFlag  : If number of minimum observations is under the 
        %                 threshold after editing process, minObsFlag takes 
        %                 value "1" else value "0". If minObsFlag takes the
        %                 value "0", updSVobj and updEDobj variables are
        %                 set to empty "[]"
        %-----------------------------------------------------------------      
                                                     
        % Outlier Detection    
          [updSVobj,updEDobj,outlierFlag,minObsFlag]=this.removeOutliers( ...
                                                     predSVobj,currEDobj);

          % Control the minumum number of observation
            if minObsFlag==1
               updSVobj=[];updEDobj=[];return;

            end
            
        
        end
       
        function [updSVobj,updEDobj,outlierFlag,minObsFlag]=removeOutliers( ...
                                                   this,SVobj,EDobj)
        %-----------------------------------------------------------------             
        % FUNCTION
        %   Function search for a possible outlier in observation set
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables  
        % OUTPUTS
        %   updEDobj : instance of 'OrbEpochData' class including updated 
        %              observations, GPS ephemerides and other auxiliary data
        %   updSVobj : instance of 'OrbStateVector' class including updated 
        %              state vector variables    
        %   outlierFlag : If any outlier detected, outlierFlag takes value 
        %                 '1', else value '0'
        %   minObsFlag  : If there exist an outlier, minObsFlag takes value 
        %                '1' else value '0'
        %----------------------------------------------------------------- 

        % Set parameters 
          stdObsNoise=this.filtSet.stat.std.measurementNoise;
          stdFactor=this.filtSet.dataEditing.outlierFactor;
        % Get apriory observations
          z=EDobj.obs.NavSol;

        % Compute expected uncertanity of the  measurements 
            [H]=this.getDesignMatrix (EDobj,SVobj);
            [rH,cH]=size(H); 
            Pp=SVobj.stateCov;
            H(:,7:end)=[];  Pp(:,7:end)=[]; Pp(7:end,:)=[];
            for i=1:length(z)
                sigmaRes(i,1)=sqrt(H(i,:)*Pp*H(i,:)'+ ...
                                         stdObsNoise(i)^2);  
            end
        % Get the predicted observations
          zc=[SVobj.position;SVobj.velocity];            
        % Search for the bad measurement
          outlierFlag=0; 
        % Compute the residuals
          res=z-zc;
        % Normalized residuals
          resNorm=res./sigmaRes;
        % RMS of normalized residuals
          rmsResNorm=std(resNorm);   
        % Compare to the threshold
          if rmsResNorm >= stdFactor
              outlierFlag=1;   
          end
        % Set the updated output epoch data and state vector objects
          updEDobj=OrbEpochData(EDobj.obs,EDobj.svn,EDobj.gpsPos, ...
                                EDobj.gpsClcCorr,EDobj.epochFateTime,...
                                EDobj.hFlag);               
          % Initialize output State Vector
            updSVobj=OrbStateVector(SVobj.stateVector, ...
                                    this.filtSet.enableDynModParam, ...
                                    this.filtSet.obsType,...
                                    SVobj.stateFateTime);    
            updSVobj.stateCov=SVobj.stateCov;                   

        % Control the observation count
          if outlierFlag==1
             minObsFlag=1;
          else
             minObsFlag=0;
          end             
            
        end
        
        function [H]=getDesignMatrix(this,EDobj,SVobj)
        %-----------------------------------------------------------------             
        % FUNCTION
        %   Function computes the design matrix for the Navigation solution  
        %   measurement model
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables       
        % OUTPUTS
        %   H     : Computed design matrix
        %            
        %----------------------------------------------------------------- 
        
        % Partial derivatives of measurement model with respect 
        % to state vector parameters 
             m=length(SVobj.stateVector);
             n=length(EDobj.obs.NavSol);
             H =zeros(n,m); H(1:n,1:6)=eye(6);
            
        end
        
        function [zc]=getCompObs(this,EDobj,SVobj)
        % FUNCTION
        %   Function computes the observations based on the measurement model 
        % INPUTS
        %   EDobj : instance of 'OrbEpochData' class including observations
        %           GPS ephemerides and other auxiliary data
        %   SVobj : instance of 'OrbStateVector' class including state
        %           vector variables       
        % OUTPUTS
        %   zc    : Computed observations
        %            
        %----------------------------------------------------------------- 
        
        % Computed Navigation Solution observable
          zc=[SVobj.position;SVobj.velocity];          
            
        end
   
    end
    
end

