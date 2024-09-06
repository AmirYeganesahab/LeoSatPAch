classdef SatelliteDynamic_v8 < handle
    % CLASS
    %   implements the satellite dynamical motions. Class includes
    %   functions for accelerations and state propogation  
    % WRITTEN BY
    %   Eren Erdogan, erdogan.eren@gmail.com
    %------------------------------------------------------------------
    
    properties  (Access=private, Hidden)
       % Parameters for Earth Gravitational model
         Cnm; Snm   % Stokes coefficients
         nmax       % maximum degree  of Stoke's coefficients
       % Parameters for Satellite
         A2m        % satellite area to mass ratio m^2/kg
         sat_mass   % mass of the satellite in unit of kg
         inc_id     % satellite orbital parameter representing inclination
                    % ( from 2 for low inclination satellites to 6 for polar
                    % orbits)  
       % Astodynamical constants    
         GM_Earth=3986004.415E+8  % Earth's gravity constant, [m^3/s^2],EGM2008
         GM_Sun= 1.32712438e+20   % Sun's gravity constant [m^3/s^2]; IAU 1976
         GM_Moon= 4902799059741.11% Moon's garvity constant
         w_earth=7.2921158553e-5  % Earth angular velocity, NIMA 1997
         AU=149597870000.0;       % Astronomical unit [m], IAU 1976
       % Object for time and reference system
         timeRefSysObj=TimeAndRefSystem();
       % Solar radiation pressure constants
         Psrad=4.56e-6       % Solar radiation pressure at 1 AU, [N/m^2], IERS 96
       % Time and Reference System Parameters
         TAI_UTC                     % Difference between TAI and UTC time
         UT1_UTC                     % Difference between UTC and TAI time
         Re =6378136.46              % Earth's mean radius(m)
         xpm; ypm;                   % polar motion parameters [in second]
       % Numeric Integration Object
         numIntegObj
       % Gravity object
         gravityObj
    end
    
    methods
      %*******************************************************************
      % Main Class
      %*******************************************************************
        function this=SatelliteDynamic_v8(A2m,inc_id,sat_mass, ...
                                       nmax,TAI_UTC,UT1_UTC,xpm,ypm)
        % Load gravity coefficients
          gravity_data= load('CS_EGM2008_n105_norm.mat');
          Cnm_full=gravity_data.C; % Normalized gravity coefficients (max 100)
          this.Cnm=Cnm_full(1:nmax+1,1:nmax+1);
          Snm_full=gravity_data.S; % Normalized gravity coefficients (max 100)
          this.Snm=Snm_full(1:nmax+1,1:nmax+1);
        % Set parameters
          this.A2m=A2m;     % satellite area to mass ratio m^2/kg
          this.sat_mass=sat_mass;
          this.nmax=nmax;   % Maximum number of Stokes coefficient degree,n
          this.inc_id=inc_id;
          this.xpm=xpm; this.ypm=ypm;   % polar motion parameters [in second]
          this.TAI_UTC=TAI_UTC;
          this.UT1_UTC=UT1_UTC;
        % Create objects
          % Numeric Integration Object 
            this.numIntegObj=NumericalIntegration('RK4');
          % Object for gravity computations
            this.gravityObj=GeoPotNormPines(this.Re,this.Cnm,this.Snm, ...
                                            this.GM_Earth,this.nmax);
        end
        
      %*******************************************************************
      % Time Propagation
      %******************************************************************* 
      
        function [varargout]=propagate(this,propType, stateModelType,position,...
                                       velocity,atmDragCoef,solarRadCoef, ...
                                       empAccelRTN,corelTime,MJD_TT0,dt,dh)
         propType='State&Transition' ;                    
        % FUNCTION
        %   Prppagetes satellite state and transition matrix through the
        %   time
        % INPUTS
        %   propType     : indicates the type. 
        %                     if type='State' only state propagation is done. 
        %                     if type='State&Sensitivity'  state vactor and state 
        %                     sensitivity matrix are propagated    
        %   stateModelType: indicates which state parameters will be
        %                   propageted. 
        %                     - If it takes "0", only position and velocity
        %                    are propageted, and the other parameters are
        %                    assumed constant thrpugh the integration.
        %                     - If it takes "1", position, velocity,
        %                     atmospheric drad coefficient, solar radiation
        %                     pressure coefficient and empirical
        %                     accelerations are integrated through the time
        %                     and the other parameters are  assumed constant
        %                     through the integration. 
        %                     - If it takes "2", position, velocity,
        %                     atmospheric drad coefficient, solar radiation
        %                     pressure coefficient,empirical accelerations
        %                     and markov process corelletion time are integrated 
        %                     time through the time, and the other parameters are
        %                     assumed constant through the integration.           
        %   position     : initial position vector
        %   velocity     : initial velocity vector
        %   atmDragCoef  : initial atmospheric drag coefficient
        %   solarRadCoef : initial solar radiation coefficient
        %   empAccelRTN  : initial empirical acceleration vector along  
        %                  radial, tangential and normal diractions
        %   corelTime    : initial corelation time for Gauss Markov Process             
        %   MJD_TT0      : modified julian date of initial time (Terrastrial Time)  
        %   dh           : step size for numerical integration (in second)
        %   dt           : total integration time (in second)
        % OUTPUTS
        %  
        %----------------------------------------------------------------  
        
        % Set the state vector
          state0=[position(:);velocity(:)];
        % Extend the state vector for propagation of dyanmical
        % model parameters and set the constant paarameters          
          if stateModelType==0; % Only position and velocity propagation
             constModelParam={atmDragCoef,solarRadCoef, ...
                             empAccelRTN,corelTime};
          elseif stateModelType==1;  
                 state0=[state0;atmDragCoef(:);solarRadCoef(:);empAccelRTN(:)];
                 constModelParam={corelTime(:)};
          elseif stateModelType==2;
                 state0=[state0;atmDragCoef(:);solarRadCoef(:);...
                         empAccelRTN(:);corelTime(:)];
                 constModelParam={[]};
          end
        % If computation of transition matrix is required, extend the state
        % vector to include initial transition matrix parameters
          n_st=length(state0); % length of state vector parameters that will 
                               % be propagated 
          if strcmp(propType,'State&Transition' )
             phi0=eye(n_st);
             state0=[state0;phi0(:)];
          end
        % Propagate state
          [state_p]=this.numIntegObj.RK4(0,state0,dt,dh,@this.derivatives ,...
                                    MJD_TT0,propType,stateModelType, ...
                                    constModelParam{:});
        % Set the function outputs
          state=state_p(1:n_st);                        % propagated state vector
          varargout{1}=state; 
          if strcmp(propType,'State&Transition' )
             phi=reshape(state_p(n_st+1:end),[n_st,n_st]); % propagated transition matrix
             varargout{2}=phi;
          end
           
        end        
        
      %*******************************************************************
      % Time Derivatives
      %*******************************************************************    
        function [d_state]=derivatives(this,dt,state,mjd_TT,propType, ...
                                        stateModelType,varargin)
        % FUNCTION
        %   Function calculates state derivatives
        % INPUTS
        %   dt    : time difference between the current and initial
        %           time (in second)
        %   state : state vector
        %   propType     : indicates the type. 
        %                  -if type='State' only state derivatives are computed. 
        %                  -if type='State&Sensitivity' partial derivatives of  
        %                   force function with respect to state parameters 
        %                   are computed
        %   stateModelType: indicates which state parameters will be
        %                   propageted. 
        %                    - If it takes "0", only derivatives of the 
        %                    position and velocity are computed, and the other 
        %                    parameters are assumed constant through the derivation.
        %                    - If it takes "1", derivatives of position, velocity,
        %                    atmospheric drad coefficient, solar radiation
        %                    pressure coefficient and empirical accelerations
        %                    are computed and the other parameters are  
        %                    assumed constant
        %                    - If it takes "2",  derivatives of position, velocity,
        %                    atmospheric drad coefficient, solar radiation
        %                    pressure coefficient,empirical accelerations
        %                    and markov process corelletion time are computed 
        %                    and the other parameters are assumed constant
        %   varragin: matlab variable type for input containing constant model
        %             paremeters. If stateModelType is "0", it is defined
        %             as follows;
        %                varragin={atmDragCoef,solarRadCoef,empAccelRTN,corelTime}
        %             If stateModelType is "1", it is defined as follows;
        %                varragin={corelTime}      
        %             If stateModelType is "2", it is defined as follows;
        %                varragin={[]}         
        %   mjd_TT : indicates the modified julian date of
        %             time (TT) at initial time.  
        % OUTPUTS
        %   d_state : derivative of given state space
        %             d_state=[vx ; vy ; vz ; ax ; ay ; az; ...]
        %----------------------------------------------------------------
        
        % Set the state vector parameters
          x=state(1) ; y=state(2) ; z=state(3);
          vx=state(4); vy=state(5); vz=state(6);
          if stateModelType==0; % Only position and velocity propagation
             atmDragCoef=varargin{1}; % Atmospheric drag coefficient
             solarRadCoef=varargin{2};% Solar radiation coefficient
             wR=varargin{3}(1);       % empirical acceleration (radial direction)
             wT=varargin{3}(2);       % empirical acceleration (tangential direction)
             wN=varargin{3}(3);       % empirical acceleration (normal direction)
             corelTime=varargin{4};   % Markov Process corellation time
          elseif stateModelType==1;  
             atmDragCoef=state(7); 
             solarRadCoef=state(8);
             wR=state(9); 
             wT=state(10); 
             wN=state(11); 
             corelTime=varargin{1}; 
          elseif stateModelType==2;
             atmDragCoef=state(7); 
             solarRadCoef=state(8);
             wR=state(9); 
             wT=state(10); 
             wN=state(11); 
             corelTime=state(12); 
          end
        % Compute Modified julian date of current epoch
          tn=mjd_TT+dt/86400;
        % Compute Rotation matrix from ECI to ECEF
          JD_TT=tn+2400000.5;	
          [year, month, day,hour,minute,second]=this.timeRefSysObj.jd2date(JD_TT);
          UTC=hour+minute/60+second/3600-(this.TAI_UTC+32.184)/3600;
          [E,dE]=this.timeRefSysObj.eci2ecef(year,month,day,UTC, ...
                                             this.xpm,this.ypm,this.TAI_UTC,...
                                             this.UT1_UTC,this.w_earth);  
        % Rotation matrix from Satellite orbital reference system to
        % geocentric system
          E1=this.timeRefSysObj.orb2geo(x,y,z,vx,vy,vz);
        % Total accelartion
          [ax,ay,az]=this.satTotalAccel(x,y,z,vx,vy,vz,atmDragCoef,solarRadCoef, ...
                                        wR,wT,wN,tn,E,E1);
        % Set the state derivatives
          d_state(1:3,1)=[vx;vy;vz];   % derivative of position vector
          d_state(4:6,1)=[ax;ay;az];   % derivative of velocity vector
          if stateModelType==1 || stateModelType==2 ;
             d_state(7,1)=0; % derivative of atmospheric drag coefficient
             d_state(8,1)=0; % derivative of solar radiation coefficient
             d_state(9:11,1)=-1/corelTime*[wR ;wT ;wN]; % derivatives of
                                                        % empirical accelerations
             if stateModelType==2
                d_state(12,1)=0; % derivative of markov process corelletain time
             end
          end
        % If both state vector and state transition matrix will be 
        % propagated,expand the derivative resulting that state derivative 
        % vector includes the components of sensitivit matrix including 
        % partial derivatives    
          if strcmp(propType,'State&Transition' )
             [F]=this.totalPartials(x,y,z,vx,vy,vz,atmDragCoef,solarRadCoef,...
                                    wR,wT,wN,corelTime,tn,E,E1,stateModelType);
             [rF,cF]=size(F);  
             phi=state(rF+1:end); phi=reshape(phi,[rF,rF]);
             dPhi=F*phi;

             %dPhi=F*eye(rF);
             d_state=[d_state;dPhi(:)];                   
          end
        end
        
      %*******************************************************************      
      % Accelerations acting on satellite 
      %*******************************************************************   
      
        function [ax,ay,az]=satTotalAccel (this,x,y,z,vx,vy,vz,Cdrag,Crad,...
                                           aR,aT,aN,mjd_tt,E,E1)
        % FUNCTION
        %   Total acceleration acting on satellite are computed at given
        %   state.
        % INPUTS
        %   x,y,z    : indicate satellite positions in Earth Fixed Sytem 
        %   vx,vy,vz : satellite velocity       
        %   mjd_tt   : Modified Julian Date of time   
        %   E        : rotation matrix from ECI to ECEF
        %   E1       : transformation matrix from satellite orbital (radial,
        %              along track, cross track coordinates)reference system
        %              to geocebtric system
        %   aR,aA,aC : empirical accelations 
        %   Cdrag    : atmospheric drag coefficient
        %   Crad     : solar radiation pressure coefficient
        % OUTPUTS
        %   ax,ay,az : tottal acceleration acting on satellite at given
        %   epoch
        %-------------------------------------------------------------

        % Accelleration in Earth fixed System due to Earth's  
        % oblateness and mass distribution 
          [geo_ax,geo_ay,geo_az]=this.gravityObj.firstDeriv(x,y,z);  
        % Acceleration due to centrifugal force
          [cent_ax,cent_ay,cent_az]=this.centrifugalAccel(x,y,z);
        % Acceleration due to coriolis effect
          [cori_ax,cori_ay,cori_az]=this.coriolisAccel(vx,vy,vz);
        % Acceleration due to atmosphric drag
            % Compute acceleration
              [drag_ax,drag_ay,drag_az]=dragAccel(this,x,y,z,...
               this.A2m,Cdrag,this.inc_id,[vx;vy;vz],mjd_tt,E);
        % Acceleration due to Lunar and Solar bodies    
            % Coordinates of Sun and Moon
              [x_sun, y_sun, z_sun]=this.sunPos(mjd_tt);
              [x_moon, y_moon, z_moon]= this.moonPos(mjd_tt);
              % Convert the Sun and Moon position to the ECEF system
                r_sun=E*[x_sun; y_sun; z_sun];
                x_sun=r_sun(1);y_sun=r_sun(2);z_sun=r_sun(3);
                r_moon=E*[x_moon; y_moon; z_moon];
                x_moon=r_moon(1);y_moon=r_moon(2);z_moon=r_moon(3);
            % Compute acceleration            
              [ls_ax, ls_ay, ls_az]=luniSolarAccel (this,...
                                x,y,z,x_sun,y_sun,z_sun, ...
                                x_moon,y_moon,z_moon);
        % Accelaration due to Solar radiation pressure
            % Compute acceleration
              [srad_ax,srad_ay,srad_az]=this.solarRadAccel(x,y,z, ...
                                      x_sun,y_sun,z_sun,Crad,this.A2m,this.Psrad);
        % Acceleration caused by solid Earth tides due to the Sun and the Moon
          [tide_ax,tide_ay,tide_az]=this.solidEartTideAccel(x,y,z, ...
              x_sun,y_sun,z_sun,x_moon,y_moon,z_moon);
        % Empirical acceleration
          [emp_accel]=E1*[aR;aT;aN];
           emp_x=emp_accel(1); emp_y=emp_accel(2); emp_z=emp_accel(3);
        % Total acceleration
          ax=geo_ax+cori_ax+cent_ax+ls_ax+srad_ax+drag_ax+tide_ax+emp_x;
          ay=geo_ay+cori_ay+cent_ay+ls_ay+srad_ay+drag_ay+tide_ay+emp_y;
          az=geo_az+cori_az+cent_az+ls_az+srad_az+drag_az+tide_az+emp_z;
            
        end              
        function [ax,ay,az]=dragAccel(this,x,y,z,A2m, ...
                                             Cdrag,inc_id,vel_rel,mjd_tt,E)
        % FUNCTION
        %   Function dragaccel computes the acceleration caused by the drag force of
        %   the atmosphere. The function computed the density of the atmosphere by
        %   Harris-Priester Density model. 
        % INPUTS 
        %   x,y,z     : position of the satellite in geocentric system
        %   Cdrag        : athmospheric drag coefficient
        %   inc_id    : parameter indicating orbit inclination        
        %   mjd       : modified julian date of time
        %   A2m       : area to mass ratio of satellite
        %   vel_rel   : relative velocity of satellite with respect to
        %               Earth' atmosphere
        %                 for example: assuming the atmosphere rotates
        %                 together with the Earth, then the relative
        %                 velocity in inertial frame;
        %                   vel_rel = v - cross(w,r)
        %                 in earth fixed frame
        %                   vel_rel = v 
        %   mjd_tt    : Modified Julian Date of terrastrial time      
        % OUTPUTS
        %   ax,ay,az  : acceleration in x,y and z direction respectively 
        %               
        %----------------------------------------------------------------
        % Atmospheric density 
          rho=this.densityHP(x,y,z,mjd_tt,inc_id,E);
        % Acceleration due to atmospheric drag
          accel=-0.5*Cdrag*A2m*rho*norm(vel_rel)*vel_rel;
          ax=accel(1); ay=accel(2); az=accel(3); 
        end
        function [rho]=densityHP(this,x,y,z,mjd,inc_id,E)
        % FUNCTION
        %   Function computes the atmospheric density using Harris Priester
        %   density model.
        % INPUTS
        %   x,y,x : satellite position in space fixed system
        %   MJD_Terrastrial: Modified Julian Date of terrastrial time  
        %   inc_id: parameter indicating orbit inclination        
        % OUTPUTS
        %   rho : atmospheric density in units of kg/m3 at given position
        % REFERANCE
        %   Satellite Orbits, Montenbruck&Gill   
        %------------------------------------------------------------
        % Set position vector
          pos=[x y z]';
        % Density coefficients
          h=[100 120 130 140 150 160 170 180 190 200 210 220 230 ...
             240 250 260 270 280 290 300 320 340 360 380 400 420 ...
             440 460 480 500 520 540 560 580 600 620 640 660 680 ...
             700 720 740 760 780 800 840 880 920 960 1000 ]';
          rhom = [497400  24900     8377      3899     2122     1263   ...
                  800.8   528.3     361.7     255.7    183.9    134.1  ...
                  99.49   74.88     57.09     44.03    34.3     26.97  ...
                  21.39   17.08     10.99     7.214    4.824    3.274  ...
                  2.249   1.558     1.091     0.7701   0.5474   0.3916 ...
                  0.2819  0.2042    0.1488    0.1092   0.0807   0.06012...
                  0.04519 0.0343    0.02632   0.02043  0.01607  0.01281...
                  0.01036 0.008496  0.007069  0.00468  0.0032   0.00221...
                  0.00156 0.00115 ]';
          rhoM = [497400   244900   8710      4059     2215     1344   ...
                  875.8    601      429.7     316.2    239.6    185.3  ...
                  145.5    115.7    93.08     75.55    61.82    50.95  ...
                  42.26    35.26    25.11     18.19    13.37    9.955  ...
                  7.492    5.684    4.355     3.362    2.612    2.042  ...
                  1.605    1.267    1.005     0.7997   0.639    0.5123 ...
                  0.4121   0.3325   0.2691    0.2185   0.1779   0.1452 ...
                  0.119    0.09776  0.08059   0.05741  0.0421   0.0313 ...
                  0.0236   0.018]';

        % Height from the reference ellipsoid
          hi = this.height(x,y,z);          
          i=0;
        % Reference Height
          while ( hi-h(i+1) > 0)
            i=i+1;
          end
        % Sun Coordinates
            [sx sy sz]=sunPos(this,mjd);
            ra_Sun  = atan2( sy, sx );                  % sun right ascension 
            dec_Sun = atan2( sz, sqrt( sx^2+sy^2));% sun declination
        % Unit vector to the apex of diurnal bulge
            lag_angle=0.523599;       %lag angle in longitude
            eb=[cos(dec_Sun)*cos(ra_Sun+lag_angle)
                cos(dec_Sun)*sin(ra_Sun+lag_angle)
                sin(dec_Sun)];        % in space fixed geocentric system
            
        % Scale Height
            Hm=(h(i)-h(i+1))/log(rhom(i+1)/rhom(i));
            HM=(h(i)-h(i+1))/log(rhoM(i+1)/rhoM(i));
        % Density at apex and antapex
            rhomh = rhom(i)*exp((h(i)-hi)/Hm);
            rhoMh = rhoM(i)*exp((h(i)-hi)/HM);
        % Diurnal density variation from apex to antapex due to solar variation
            er=E'*pos;                  % in space fixed geocentric system        
            er = (er/norm(er));
            rho=rhomh+(rhoMh-rhomh)*((1/2+ dot(er,eb)/2)^(inc_id/2));
            rho=rho/(10^12); %density from gr/km3 to kg/m3
        end   
        
        function [h] = height(this,x,y,z)
        % FUNCTION
        %   Determines the ellipsoidal height
        % INPUTS
        %   x,y,z : position
        % OUTPUTS
        %   h=ellipsoidal height in km
        %---------------------------------------------------------

        a=6378137;                 % WGS84 semimajor axis
        e2=0.00669438;             % First eccentricity
        e_2=0.006739496;           % Second eccentricity
        b=sqrt(a^2/(1+e2));        % WGS84 semiminor axis

        pos=[x y z]';
        L=atan2(pos(2),pos(1));           % elipsoidal longitude
        p=sqrt(pos(1)^2+pos(2)^2);        
        nu=atan2(pos(3)*a,p*b);
        B=atan2((pos(3)+e_2*(sin(nu))^3),...
                (p-e2*a*(cos(nu)^3)));    % elipsoidal latitude
        h=(p*cos(B)+pos(3)*sin(B)-a ...
           *sqrt(1-e2*(sin(B))^2))/1000;  % elipsoidal height in km
        end      
        
        function [ax ay az]=luniSolarAccel(this,x,y,z, ...
                            x_sun,y_sun,z_sun,x_moon,y_moon,z_moon)
        %*******************************************************
        %FUNCTION
        %   Function luniSolarAccel returns the low precision acceleration caused
        %   by sun and moon 
        %INPUTS
        %   x,y,z : satellite receiver Position in space fixed system
        %OUTPUTS
        %   ax ay az: total acceleration due to mass of sun and moon
        %        e.g. msAccel[ax,ay,az]
        %REFERANCE
        %   Montenbruck and Gill, Satellite Orbits
        %--------------------------------------------------------------

        % acceleration of satellite due to Sun
            [ax_sun ay_sun az_sun]=this.pointMassAccel(x,y,z, ...
                                   x_sun,y_sun,z_sun,this.GM_Sun);
        % acceleration of satellite due to Moon
            [ax_moon ay_moon az_moon]=this.pointMassAccel(x,y,z, ...
                                      x_moon,y_moon,z_moon,this.GM_Moon);
        % total acceleration by Sun and Moon
            ax=ax_sun+ax_moon;
            ay=ay_sun+ay_moon;
            az=az_sun+az_moon;

        end   
        
        function [ax ay az]=solarRadAccel (this,x,y,z,x_sun,y_sun,z_sun,Crad,A2m,P )
        %*******************************************************
        %FUNCTION
        %   solarRadAccel computes the acceleration acting on satellite 
        %   due to direct solar radiation. satellite is assumed symmetrically
        %   perfect and rotationally invariant sphere with uniform optical properties 
        %   so called cannonball approach. Cylindrical shadow model is used
        %   to find whether satellite is in Earth's shadow.
        %INPUTS
        %   x,y,z: satellite coordinates in geocentric system
        %   x_sun,y_sun,z_sun : sun position in geocentric reference system
        %   Crad: radiation pressure coefficient
        %   A2m       : area to mass ratio of satellite
        %   P         : solar radiation pressure at 1AU
        %OUTPUTS
        %   ax ay az: total acceleration 
        %REFERANCE
        %--------------------------------------------------------------
        % Constract the vectors
          r_sat=[x;y;z];             % satellite position vector 
          r_sun=[x_sun;y_sun;z_sun]; % Sun position vector
        % Compute shadow parameter
          v=this.cylindricalShadowModel(r_sat,r_sun);
        % Relative position vector of Sun
          n=r_sat-r_sun;
        % total acceleration
          accel=v*P*Crad*A2m*this.AU^2/(norm(n)^3)*n;
          ax=accel(1);
          ay=accel(2);
          az=accel(3);
        end 
        
        function [v]=cylindricalShadowModel(this,r_sat,r_sun)
        %*******************************************************
        %FUNCTION
        %   Computes the fact that whether satellte is in sunlight or 
        %   Erath's shadow 
        %INPUTS
        %   r_sat : satellite position vector in geocentric system
        %   r_sun : sun position vector in geocentric reference system
        %OUTPUTS
        %   v : Shadow parameter, v, takes the value one if satellite is in
        %       sunlight else it takes value zero 
        %--------------------------------------------------------------
        
        % Compute the cases
          s1=r_sat'*r_sun; % s1=cos(theta)
          s2 = norm(cross(r_sat,r_sun))/norm(r_sun); % s2=sin(theta)
        % Find the shadow parameter, v
          if s1>0 || s2>this.Re; v=1; else v=0; end
        end
        
        function [ax,ay,az]=coriolisAccel(this,vx,vy,vz)
        % FUNCTION
        %   Determines acceleration due to coriolis effect 
        % INPUTS
        %   vx,vy,vz  : satellite velocity
        % OUTPUTS
        %   ax,ay,az  : acceleraion due to coriolis effet
        %-----------------------------------------------------------------
          % instantaneous polar motion matrix
            pol_m=this.timeRefSysObj.rot('y',this.xpm/3600)* ...
                  this.timeRefSysObj.rot('x',this.ypm/3600);  
          % velocity vector
            vel=[vx,vy,vz]';
          % Approximate Earth rotation vector
            omega_earth=pol_m*[0.0; 0.0; this.w_earth];
          % Coriolis acceleration
            accel=-2*cross(omega_earth,vel);
            ax=accel(1); ay=accel(2); az=accel(3);
        end
        
        function [ax,ay,az]=centrifugalAccel(this,x,y,z)
        % FUNCTION
        %   Determines acceleration due to centrifigal effect considered in
        %   non-inertial reference system
        % INPUTS
        %   x,y,z  : satellite position
        % OUTPUTS
        %   ax,ay,az  : acceleraion due to centrifigal effet
        %-----------------------------------------------------------------
          % instantaneous polar motion matrix
            pol_m=this.timeRefSysObj.rot('y',this.xpm/3600)* ...
                  this.timeRefSysObj.rot('x',this.ypm/3600);  
          % position vector
            pos=[x y z]';
          % Approximate Earth rotation vector
            omega_earth=pol_m*[0.0; 0.0; this.w_earth];
          % Centrifugal Acceleration
            accel=-cross(omega_earth,(cross(omega_earth,pos)));
            ax=accel(1); ay=accel(2); az=accel(3);

        end
        
        function [ax,ay,az]=solidEartTideAccel(this,x,y,z,x_sun,y_sun,z_sun, ...
                                               x_moon,y_moon,z_moon)
        % FUNCTION
        %   Function determines accelerations of the satellite caused by 
        %   solid Earth tides due to the Sun and the Moon
        % INPUTS
        %   x,y,z                : satellite position
        %   x_sun,y_sun,z_sun    : coordinates of the Sun
        %   x_moon,y_moon,z_moon : coordinates of the Moon
        %REFERANCE
        %   Satellite Geodesy,2nd Editon, Günter Seeber, p.101
        %   Force Modelling of GPS satellite orbits,Rizos C.,Stolz A. 1985       
        % OUTPUTS
        %   ax,ay,az  : the total acceleration
        %----------------------------------------------------------------- 
        % Set constants
          k2=0.3;  % approximate second degree love number
        % Set coordinate vectros
          r_sat=[x y z]'; r_sat_norm=norm(r_sat);
          r_sun=[x_sun,y_sun,z_sun]'; r_sun_norm=norm(r_sun);
          r_moon=[x_moon,y_moon,z_moon]'; r_moon_norm=norm(r_moon);
        % Cosinus of angels between position vectors of satellite and planetary body
          cos_sat_sun=(r_sat'*r_sun)/(r_sat_norm*r_sun_norm);
          cos_sat_moon=(r_sat'*r_moon)/(r_sat_norm*r_moon_norm);
        % Tidal acceleration due to sun and moon
          fs=k2/2*(this.Re^5/r_sat_norm^4)*(this.GM_Sun/r_sun_norm^3);
          a_sun=fs*( (3-15*cos_sat_sun^2)*r_sat/r_sat_norm+...
                (6*cos_sat_sun)*(r_sun/r_sun_norm) );
          fm=k2/2*(this.Re^5/r_sat_norm^4)*(this.GM_Moon/r_moon_norm^3);
          a_moon=fm*( (3-15*cos_sat_moon^2)*r_sat/r_sat_norm+...
                (6*cos_sat_moon)*(r_moon/r_moon_norm) );    
          a=a_sun+a_moon;
        % Set outputs
          ax=a(1);ay=a(2);az=a(3);
        end
        
        function [ax ay az]=pointMassAccel(this,x_sat,y_sat,z_sat,...
                            a_mass,y_mass,z_mass,gravCoeff)
        %*******************************************************
        %PURPOSE
        %   Function ppointMassAccel returns the gravitational accelarion of the 
        %   celestial mass on satellite. It is assumed that celestial body is in 
        %   the form of point mass.
        %REFERANCE
        %   Satellite Orbits, Montenbruck&Gill, p. 69
        %INPUTS
        %   x_sat,y_sat,z_sat : Satellite coordinatess, [x,y,z] 
        %   a_mass,y_mass,z_mass : Celestial body coordinates, [x,y,z] 
        %   gravCoeff : Gravitational coefficient of point mass (GM)
        %OUTPUTS
        %   accel: includes components of accelaration on satellite due to point
        %   mass.
        %        e.g. accel[ax,ay,az]
        %********************************************************
        coordMass=[a_mass y_mass z_mass]';
        coordSat=[x_sat y_sat z_sat]';
        sr=coordMass-coordSat;
        accel=gravCoeff*( sr/(power(norm(sr),3))- coordMass/(power(norm(coordMass),3)));
        ax=accel(1); ay=accel(2); az=accel(3);

        end     
                
        
      %*******************************************************************        
      % Partial derivatives 
      %*******************************************************************      
        function [F]=totalPartials(this,x,y,z,vx,vy,vz,Cdrag,Csrad,...
                                   wR,wT,wN,cT,mjd_tt,E,E1,stateModelType)
        % FUNCTION
        %   Computes total partial derivatives of equation of motion with
        %   respect to state vector components
        % INPUTS
        %   x,y,z    : satellite positions in erath fixed sytem
        %   vx,vy,vz : satellite velocity in erath fixed sytem 
        %   wR,wT,wN : empirical accelerations in radial, tangential and normal 
        %              directions, respectively
        %   Cdrag    : atmospheric drag coefficient
        %   Csrad    : solar radiation coefficient     
        %   cT       : markov process corellation time
        %   E        : rotation matrix from ECI to ECEF
        %   E1       : transformation matrix from satellite orbital (radial,
        %              along track, cross track coordinates)reference system
        %              to geocentric system
        %   mjd_tt   : Modified julian date of terrastrial time
        %   stateModelType:  - If it takes "0", only partial derivatives of  
        %                    the position and velocity are computed.
        %                    - If it takes "1", partial derivatives of position, 
        %                    velocity, atmospheric drad coefficient, solar 
        %                    radiation pressure coefficient and empirical 
        %                    accelerations are computed.
        %                    - If it takes "2",  partial derivatives of position, 
        %                    velocity, atmospheric drad coefficient, solar 
        %                    radiation pressure coefficient,empirical accelerations
        %                    and markov process corelletion time are computed 
        % OUTPUTS
        %   total_partials : matrix including partial derivatives
        %--------------------------------------------------------------
        
        % Partial Derivatives of gravity acceleration 
          [geo_da_dr]=this.gravityObj.secondDeriv(x,y,z);
        % Partial Derivatives of coriolis acceleration   
          [cori_da_dv]=this.coriolisAccelPartials();
        % Partial Derivatives of centrifugal acceleration           
          [cent_da_dr]=this.centrifugalAccelPartials();
        % Partial derivatives of solar radiation pressure acceleration
          % Coordinates of Sun
             [x_sun, y_sun, z_sun]= this.sunPos(mjd_tt);  
              % Convert the Sun and Moon position to the ECEF system
                r_sun=E*[x_sun; y_sun; z_sun];
                x_sun=r_sun(1);y_sun=r_sun(2);z_sun=r_sun(3);
          % Partial derivatives
             [srad_da_dr,da_dCsrad]=this.solarRadAccelPartials(...
              x,y,z,x_sun,y_sun,z_sun,Csrad,this.A2m,this.Psrad);    
        % Partial Derivatives of drag acceleration 
           [drag_da_dr,drag_da_dv,da_dCdrag]=dragAccelPartialsEF(this,mjd_tt,...
           x,y,z,vx,vy,vz,Cdrag,this.inc_id,this.A2m,E);
        % Partial Derivatives for the empirical acceleration components
           [emp_da_dr,emp_da_dv,da_dw,dWdot_dcT,dWdot_dw]=this.empAccelPartials(...
                                                x,y,z,vx,vy,vz,wR,wT,wN,cT,E1);         
        % Construct sensitivity matrix
          % Create the sensitivity matrix 
            F=zeros(6); 
            F(1:3,4:6)=eye(3);                                      % dv/dv  
            F(4:6,1:3)=geo_da_dr+cent_da_dr+drag_da_dr+srad_da_dr;  % da/dr
            F(4:6,4:6)=cori_da_dv+drag_da_dv;                       % da/dv
            if stateModelType==1 || stateModelType==2
               F(4:6,7)=da_dCdrag; F(4:6,8)=da_dCsrad; 
               F(4:6,9:11)=da_dw; F(9:11,9:11)=dWdot_dw;
               F(4:6,1:3)=F(4:6,1:3)+emp_da_dr;  
               F(4:6,4:6)=F(4:6,4:6)+emp_da_dv;
               if stateModelType==2
                  F(9:11,12)=dWdot_dcT;
                  F(12,:)=0;
               end
            end          
        end
        
        function [da_dr]=geoAccelPartials(this,x,y,z,Re,GM,C,S,nmax,mmax)
        % FUNCTION
        %   Computes the partial derivatives of gravity acceleration with
        %   respect to state matrix. Gravity acceleration is not dependent on
        %   velocity. Thus partials with respect to velocity component is zero
        % INPUTS
        %   x,y,z : satellite positions in space fixed system  
        %   Re    : Earth radius
        %   GM    : gravitational constant of Earth
        %   C,S   : stoke's coefficients
        %   nmax,mmax  : maximum degree of stoke's coefficients         
        % OUTPUTS
        %   geo_partials:  matrix including partial derivatives
        %                 [dax/dx   dax/dy   dax/dz  
        %   geo_partials=  day/dx   day/dy   day/dz 
        %                  daz/dx   daz/dy   daz/dz  ]        
        %
        %------------------------------------------------------------

        % Partials with respect to position
          dax_dx =0;   dax_dy=0; dax_dz=0; day_dz=0; daz_dz=0;
          [V,W]=this.VW(x,y,z,nmax+2,nmax+2,Re);
          for m=0:mmax
              mm=m+1;
              for n=m:nmax
                  nn=n+1;
                  if m==0
                     dax_dx = dax_dx+0.5*((C(nn,1)*V(nn+2,3))+ ...
                              (n+1)*(n+2)*(C(nn,1)*V(nn+2,1)));
                     dax_dy = dax_dy+0.5*(C(nn,1)*W(nn+2,3));
                     dax_dz = dax_dz+(n+1)*(C(nn,1)*V(nn+2,2));
                     day_dz = day_dz+(n+1)*(C(nn,1)*W(nn+2,2));
                  end
                  if m==1
                     dax_dx = dax_dx+0.25*((C(nn,2)*V(nn+2,4) ...
                              +S(nn,2)*W(nn+2,4))+n*(n+1)*(-3*C(nn,2) ...
                               *V(nn+2,2)-S(nn,2)*W(nn+2,2)));
                     dax_dy = dax_dy+0.25*((C(nn,2)*W(nn+2,4)- ...
                              S(nn,2)*W(nn+2,4))+n*(n+1)*(-C(nn,2)* ...
                              W(nn+2,2)-S(nn,2)*V(nn+2,2)));
                  end
                  if m>1
                     dax_dx = dax_dx+ ...
                              0.25*((C(nn,mm)*V(nn+2,mm+2)+S(nn,mm)...
                              *W(nn+2,mm+2))+2*(n-m+1)*(n-m+2)*(-C(nn,mm) ...
                              *V(nn+2,mm)-S(nn,mm)*W(nn+2,mm)) +(n-m+4)...
                              *(n-m+3)*(n-m+2)*(n-m+1)*(C(nn,mm)*V(nn+2,mm-2)...
                              +S(nn,mm)*W(nn+2,mm-2)));
                     dax_dy = dax_dy+0.25*((C(nn,mm)*W(nn+2,mm+2)...
                              -S(nn,mm)*V(nn+2,mm+2))+(n-m+4)*(n-m+3) ...
                              *(n-m+2)*(n-m+1)*(-C(nn,mm)*W(nn+2,mm-2)...
                              +S(nn,mm)*W(nn+2,mm-2)));
                  end
                  if m>0
                     dax_dz = dax_dz+((n-m+1)/2*(C(nn,mm)...
                              *V(nn+2,mm+1)+S(nn,mm)*W(nn+2,mm+1))...
                              +(n-m+3)*(n-m+2)*(n-m+1)/2*(-C(nn,mm) ...
                              *V(nn+2,mm-1)-S(nn,mm)*W(nn+2,mm-1)));
                     day_dz = day_dz+((n-m+1)/2*(C(nn,mm)...
                              *W(nn+2,mm+1)-S(nn,mm)*V(nn+2,mm+1))...
                              +(n-m+3)*(n-m+2)*(n-m+1)/2*(C(nn,mm)...
                              *W(nn+2,mm-1)-S(nn,mm)*V(nn+2,mm-1)));
                  end
                  daz_dz = daz_dz+((n-m+2)*(n-m+1) ...
                           *(C(nn,mm)*V(nn+2,mm)+S(nn,mm)*W(nn+2,mm)));
               end
          end
        % Sum of diagonal elements of partial matrix is zero
            day_dy = -dax_dx-daz_dz;
        % matrix including partials is symmetrical
            day_dx = dax_dy;
            daz_dx = dax_dz;
            daz_dy = day_dz;
        % Partials with respect to position
            da_dr=(GM/(Re^3)).*[dax_dx  dax_dy  dax_dz
                                day_dx  day_dy  day_dz
                                daz_dx  daz_dy  daz_dz];

        end    
        
        function [da_dr,da_dv,da_dCdrag]=dragAccelPartialsEF(this,mjd_tt,...
                                       x,y,z,vx,vy,vz,Cdrag,inc_id,A2m,E)
        % FUNCTION
        %   Computes the partial derivatives of drag acceleration with
        %   respect to position velocity vectors of satellite and drag 
        %   coefficient. Computation is carried out in Earth fixed system.
        % INPUTS
        %   x,y,z    : satellite position in Earth fixed system        
        %   vx,vy,vz : satellite velocitiy in Earth fixed system 
        %   Cdrag       : atmospheric drag coefficient
        %   inc_id   : parameter indicating orbit inclination
        %   A2m      : satellite area to mass ration in kg/m^2
        %   mjd_tt : Modified julian date of terrastrial time            
        % OUTPUTS
        %   da_dr  : partials with respect to position 
        %   da_dv  : partials with respect to velocity 
        %   da_dCdrag : partials with respevt to deag coefficients
        %----------------------------------------------------------------
          vel=[vx vy vz]';
          vel_norm=norm(vel);
        % Atmospheric density 
          rho=densityHP(this,x,y,z,mjd_tt,inc_id,E);
        % Compute da/dv
          da_dv=-0.5*Cdrag*A2m*rho*((vel*vel')/vel_norm+vel_norm*eye(3)); 
        % Compute dp/dr by different quotient approximation
          x1=x+1000; y1=y+1000; z1=z+1000;
          rho1=densityHP(this,x1,y1,z1,mjd_tt,inc_id,E);
          dp_dr=(rho1-rho)./[1000,1000,1000];
        % Compute da/dr
          da_dr=-0.5*Cdrag*A2m*vel_norm*vel*dp_dr;
        % Compute Sensitivity, da/dCD
          da_dCdrag=-0.5*A2m*rho*vel_norm*vel;
          
        end   
        
        function [da_dr,da_dCsrad]=solarRadAccelPartials(this,x,y,z,...
                                            x_sun,y_sun,z_sun,Csrad,A2m,P)
        %FUNCTION
        %   Computes the partial derivatives of acceleration due to solar
        %   radiation pressure with respect to satellite position and
        %   radiation pressure coefficient
        %           da_dr , da_dCrad
        %INPUTS
        %   x,y,z: satellite coordinates in geocentric system
        %   x_sun,y_sun,z_sun : sun position in geocentric reference system
        %   Crad: radiation pressure coefficient
        %   A2m       : area to mass ratio of satellite
        %   P         : solar radiation pressure at 1AU
        %OUTPUTS
        %   ax ay az: total acceleration 
        %REFERANCE
        %--------------------------------------------------------------
        
        % Compute auxiliary components
          r_sat=[x;y;z];             % satellite position vector 
          r_sun=[x_sun;y_sun;z_sun]; % Sun position vector
          n=r_sat-r_sun; n_norm=norm(n); n_norm3=n_norm^3;n_norm5=n_norm^5;
          k=1/n_norm3*eye(3)-3*n*(n'/n_norm5);
        % partial derivatives with respect to satellite position vector
          da_dr=P*Csrad*A2m*this.AU^2*k;
        % partial derivatives with respect to radiation pressure
        % coefficient
          da_dCsrad= P*A2m*this.AU^2*n./n_norm3;   
        end
        
        function [da_dv]=coriolisAccelPartials(this)
        % FUNCTION
        %   Computes the partial derivatives of coriolis acceleration with
        %   respect to state matrix. Coriaolis is depend on only velocity 
        %   of the moving object in a rotating frame.   
        % INPUTS
        % OUTPUTS
        %   da_dv    : partial derivatives of coriolis acceleration with
        %              respect to velocity
        %-------------------------------------------------------------
        % instantaneous polar motion matrix
          pol_m=this.timeRefSysObj.rot('y',this.xpm/3600)* ...
                this.timeRefSysObj.rot('x',this.ypm/3600);  
        % Approximate Earth rotation vector
          omega_earth=pol_m*[0.0; 0.0; this.w_earth];
          wx=omega_earth(1); wy=omega_earth(2); wz=omega_earth(3);
        % Partial derivatives with respect to velocity
          Rw= [0    -wz    wy
               wz    0    -wx
              -wy    wx    0];
          da_dv=-2*Rw';
        end
        
        function [da_dr]=centrifugalAccelPartials(this)
        % FUNCTION
        %   Computes the partial derivatives of centrifugal acceleration with
        %   respect to state vector. Centrifugal is depend on only position 
        %   of the moving object in a rotating frame.   
        % INPUTS
        % OUTPUTS
        %   da_dr    : partial derivatives of centrifugal acceleration with
        %              respect to position
        %-------------------------------------------------------------
        % instantaneous polar motion matrix
          pol_m=this.timeRefSysObj.rot('y',this.xpm/3600)* ...
                this.timeRefSysObj.rot('x',this.ypm/3600);  
        % Approximate Earth rotation vector
          omega_earth=pol_m*[0.0; 0.0; this.w_earth];
          wx=omega_earth(1); wy=omega_earth(2); wz=omega_earth(3);
        % Partial derivatives with respect to velocity
          Rw= [0    -wz    wy  % cross-product matrix
               wz    0    -wx
              -wy    wx    0];
          da_dr=-(Rw*Rw)';
        end        
        
        function [da_dr,da_dv,da_w,dWdot_dcT,dWdot_dw]=empAccelPartials(this,x,y,z, ...
                                                       vx,vy,vz,wR,wT,wN,cT,E)
        % FUNCTION
        %   Computes the partial derivatives of empirical accelerations with
        %   respect to state vector.   
        % INPUTS
        %   x,y,z    : satellite positions in erath fixed sytem
        %   vx,vy,vz : satellite velocity in erath fixed sytem         
        %   wR,wT,wN : empirical accelerations in radial, tangential and normal 
        %              directions, respectively
        %   cT       : Markov process corelletaion time
        %   E        : rotation matrix from orbital system to Earth fixed system
        % OUTPUTS
        %   da_dr    : partial derivatives of acceleration with
        %              respect to position vector
        %   da_dv    : partial derivatives of acceleration with
        %              respect to velocity vector
        %   dWdot_dw : partial derivatives of  acceleration with
        %              respect to empirical acceleration component  
        %   dWdot_dcT: partial derivatives of  emprical acceleration with
        %              respect to markov process corellation coefficient
        %-------------------------------------------------------------            
         
        % Compute the unit vector partial derivatives of the rotation matrix
        % E=(eR,eT,eN)
            eR=E(:,1); eT=E(:,2); eN=E(:,3); % unit vector
            r=[x y z]'; v=[vx vy vz]';       % position and velocity vector
            r_norm=norm(r); rxv_norm=norm(cross(r,v));
            Xr= [0 -z y ; z 0 -x ; -y x 0];
            Xv= [0 -vz vy ; vz 0 -vx ; -vy vx 0];
            XeR= [0 -eR(3) eR(2) ; eR(3) 0 -eR(1) ; -eR(2) eR(1) 0];
            XeN= [0 -eN(3) eN(2) ; eN(3) 0 -eN(1) ; -eN(2) eN(1) 0];
          % Partial derivatives of unit vectors with respect to position vector
            deR_dr=-1/r_norm*(eye(3)+eR*eR');
            deN_dr=-1/rxv_norm*(eye(3)-eN*eN')*Xv;
            deT_dr=-XeR*deN_dr+XeN*deR_dr;   
          % Partial derivatives of unit vectors with respect to velocity vector
            deR_dv=-zeros(3);
            deN_dv=-1/rxv_norm*(eye(3)-eN*eN')*Xr;
            deT_dv=-XeR*deN_dv;
        % Total partial derivatives of accelerations with respect to
        % position vector
          da_dr=wR*deR_dr + wT*deT_dr + wN*deN_dr;
        % Total partial derivatives of accelerations with respect to
        % velocity vector
          da_dv=wR*deR_dv + wT*deT_dv + wN*deN_dv;
        % Partial derivatives with respect to empirical acceleration
          da_w=E;
        % Partial derivatives of  emprical acceleration with respect to 
        %  markov process corellation time coefficient         
          dWdot_dcT=-[wR;wT;wN];
        % Partial derivatives of  emprical accelerations with respect to 
        % its components         
          dWdot_dw=-1/cT*eye(3);
            
        end
        
      %*******************************************************************        
      % Position of Sun and Moon
      %*******************************************************************
      
        function [x y z]= sunPos(this,MJD_Terrastrial)
        %*******************************************************
        % FUNCTION
        %   Function sunPos returns the low precision Solar coordinates,determined
        %   via simplied analytical formula, with respect to equinox and ecliptic of 
        %   epoch regarding to J2000 
        % INPUTS
        %   MJD_Terrastrial: Modified Julian Date of terrastrial time
        % OUTPUTS
        %   posSun: includes component of position vector x,y and z
        %        refering to mean equinox and ecliptic of J2000.
        %        e.g. pos[x,y,z]
        % REFERANCE
        %   Satellite Orbits, Montenbruck&Gill, p. 70-71
        %----------------------------------------------------------------
        % Modified Julian date of J2000
          MJD_J2000=2451545.0-2400000.5;       
        % julian centuries from epoch J2000
          t=(MJD_Terrastrial-MJD_J2000)/36525.0; 
        % Orbital elements 
          % Mean Anomaly in radian
            aux=0.9931267 + 99.9973583*t;
            aux=aux-floor(aux);
            M=2*pi*aux; 
          % Sun's ecliptic longitude
            aux=0.7859444+M/(2*pi)+(6892.0*sin(M)+72.0*sin(2.0*M))/ 1296.0e3;
            aux=aux-floor(aux);
            L = 2*pi* aux;  
          % Distance in meter
            r = 149.619e9 - 2.499e9*cos(M) - 0.021e9*cos(2*M);              
        % Position vector with respect to equator
          % cartesian coordinates in ecliptic system, in meter
            posSun=[r*cos(L),r*sin(L),0.0]';   
          % obliquity of the ecliptic, in decimal degree
            e_o=23.43929111;       
          % position vector with respect to equator
            posSun=this.timeRefSysObj.rot('x',-e_o)*posSun;   
            x=posSun(1); y=posSun(2);z=posSun(3);
        end    
        
        function [x y z]= moonPos(this,MJD_Terrastrial)
        % FUNCTION
        %   Function posMoon returns the low precision Lunar coordinates,determined
        %   via simplied analytical formula, with respect to equinox and ecliptic of 
        %   epoch regarding to J2000 
        % INPUTS
        %   MJD_Terrastrial: Modified Julian Date of epoch
        % OUTPUTS
        %   [x y z]: low precision position of moon
        % REFERANCE
        %   Montenbruck, Satellite Orbits
        %-----------------------------------------------------------------
        
        % MJD of J2000
          MJD_J2000=2451545.0-2400000.5;    
        % julian centuries from epoch J2000    
          t=(MJD_Terrastrial-MJD_J2000)/36525.0;              
        % Fundamental Arguments    
          % Mean longitude of the moon in radian
            aux=0.606433 + 1336.851344*t ;
            Lo =aux-floor(aux);   
          % Moon's mean anomaly in radian
            aux=0.374897 + 1325.552410*t ;
            lm = 2*pi*(aux-floor(aux)); 
          % Sun's mean anomaly in radian
            aux=0.993133 + 99.997361*t;
            ls = 2*pi*(aux-floor(aux));  
          % Mean angular distance of Moon from the ascending node F in radian
            aux=0.259086 + 1342.227825*t;
            F  = 2*pi*(aux-floor(aux));     
          % Difference between the mean longitudes of the Sun and Moon in radian
            aux=0.827361 + 1236.853086*t;
            D  = 2*pi*(aux-floor(aux));               
        % Moon's longitude,L, with respect to equinox and ecliptic of the year 2000    
            dL = +22640*sin(lm) - 4586*sin(lm-2*D) + 2370*sin(2*D) +  769*sin(2*lm) ...
               -668*sin(ls) - 412*sin(2*F) - 212*sin(2*lm-2*D) - 206*sin(lm+ls-2*D) ...
               +192*sin(lm+2*D) - 165*sin(ls-2*D) - 125*sin(D) - 110*sin(lm+ls) ...
               +148*sin(lm-ls) - 55*sin(2*F-2*D);
            aux=Lo + dL/1296.0e3;
            L=2*pi*(aux-floor(aux));              % in radian
        %Moon's latitude,B,
            S  = F + (dL+412*sin(2*F)+541*sin(ls))/(3600.0*180.0/pi);
            h  = F-2*D;
            N  = -526*sin(h) + 44*sin(lm+h) - 31*sin(-lm+h) - 23*sin(ls+h) ... 
                 +11*sin(-ls+h) - 25*sin(-2*lm+F) + 21*sin(-lm+F);
            B = ( 18520.0*sin(S) + N ) / (3600.0*180.0/pi);         % in radian

        % Moon's Distance from center of Earth,r,
            r = 385000e3-20905e3*cos(lm)-3699e3*cos(2*D-lm)...
                -2956e3*cos(2*D)-570e3*cos(2*lm)+246e3*cos(2*lm-2*D)... 
                -205e3*cos(ls-2*D)-171e3*cos(lm+2*D)...
                -152e3*cos(lm+ls-2*D); % distance in meter
        % Position vector with respect to equator
          % cartesian coordinates
            posMoon=[r*cos(L)*cos(B); r*sin(L)*cos(B);r*sin(B)];   
          % obliquity of the ecliptic
            e_o=23.43929111;     
          % position vector with respect to equator
            posMoon=this.timeRefSysObj.rot('x',-e_o)*posMoon;         
            x=posMoon(1); y=posMoon(2); z=posMoon(3);

        end

        
    end
    
end

