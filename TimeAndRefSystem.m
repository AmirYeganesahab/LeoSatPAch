classdef TimeAndRefSystem < handle
 
    
    properties
    end
    
    methods
        
        function this=TimeAndRefSystem()
        end

        function [za,thetaa,zetaa] = precession(this,t)
        % PRECESSION returns the 3 precession angles
        %
        % HOW [za,thetaa,zetaa] = precession(t)
        % IN  t      - time since epoch J2000.0 [in Julian centuries (of 36525 days)]
        % OUT za    \
        %     thetaa - precession angles        [decimal deg.]
        %     zetaa /
        %
        % NB  The accuracy of the formulae is about 1 arcsec. 

        %------------------------------------------------------------------
        % Sneeuw/Zebhauser, IAPG, TU Munich                       04/01/96
        %------------------------------------------------------------------
        % uses none
        %
        %------------------------------------------------------------------
        % revision history
        % .Oct.2001,NS: translation of original PREZWINK.M 
        %------------------------------------------------------------------
        % remarks
        %
        %------------------------------------------------------------------
        % Here we go
        zetaa  = (2306.2181 * t + 0.30188 * t.^2 + 0.017998 * t.^3)/3600;
        za     = (2306.2181 * t + 1.09468 * t.^2 + 0.18203 * t.^3)/3600;
        thetaa = (2004.3109 * t - 0.42665 * t.^2 - 0.041833 * t.^3)/3600;  
        end
        function [eps0,deps,dpsi] = nutation(this,t)
        % NUTATION returns the 3 nutation angles
        %
        % HOW [eps0,deps,dpsi] = nutation(t)
        % IN  t    - time since epoch J2000.0   [in Julian centuries (of 36525 days)]
        % OUT eps0 - mean obliquity of ecliptic [decimal deg.]
        %     deps - nutation in obliquity      [decimal deg.]
        %     dpsi - nutation in longitude      [decimal deg.]
        %
        % NB  The accuracy of the formulae used is about 1 arcsec.

        %------------------------------------------------------------------
        % Sneeuw/Zebhauser, IAPG, TU Munich                        04/01/96
        %------------------------------------------------------------------
        % uses none
        %------------------------------------------------------------------
        % revision history
        % .Oct.2001,NS: translation of original NUTWINK.M 
        %------------------------------------------------------------------
        % remarks
        %------------------------------------------------------------------
        
        % Time difference since J2000
        d  =  t * 36525;
        % Auxiliary functions
        f1 = (125 - 0.05295 * d) / 180 * pi;
        f2 = (200.9 + 1.97129 * d) / 180 * pi;
        % Here we go
        eps0   = (84381.448 - 46.815 * t)/3600;
        dpsi = (-0.0048 * sin(f1) - 0.0004 * sin(f2));
        deps = ( 0.0026 * cos(f1) + 0.0002 * cos(f2));
        end    
        function [R,gast]= earthRotation (this,t,ut1,type)
        %--------------------------------------------------------
        % gast : Greenwich apparent sideral time 
        % R : % Rotation matrix for gast
        %
        % t: time since epoch J2000.0, in form of Julian centuries
        % ut1 : universal time (ut1) at epoch, in hour
        % type : takes the value 'gast' or 'gmst'
        %--------------------------------------------------------

        % Rotating position vector with gast
            [gast,gmst]=jul2gast(t,ut1);
        % Rotation matrix for gast
            if strcmp(type,'gast')
                R=this.rot('z',-gast);
            elseif strcmp(type,'gmst')
                R=this.rot('z',-gmst);
            else
                error ('Wrong type. It takes gast or gmst')
            end
        end     
        function [E,dE]=eci2ecef(this,y,m,d,h,xp,yp,TAI_UTC,UT1_UTC,w_ecef)
        % FUNCTON
        %   Computes rotation matrix and its time derivative from 
        %   earth centered inertial frame(ECI, especially j2000 sytem) to 
        %   earth centered earth fixed reference frame (ECEF)
        %
        %   In computation of time derivative of rotation matrix,
        %   nutation,precession and time derivative assumed constant. Thus
        %     dE/dt=Pm* dEr/dt * N *P 
        %       where Pm: rotation matrix de to polar motion
        %             N : rotation matrix de to nutation
        %             P : rotation matrix de to precession
        %             Er : rotation matrix de to earth rotation
        %             
        % INPUTS
        %   y   : 4 digit year 
        %   m   : month
        %   d   : day
        %   h   : UTC (decimal hours)
        %   xp,yp  : polar motion parameters [in second]
        %   TAI_UTC : differnce between TAI and UTC [in second]
        %   UT1_UTC : differnce between UT1 and UTC [in second]
        %   w_ecef : angular velocity of ECEF system 
        % OUTPUTS
        %   E   : rotation matrix from icrf to itrf
        %   dE  : first order time derivative of rotation matrix
        % REFERENCES
        %   reference : Satellite Orbits, Montenbruck&Gill
        %   written by : Eren Erdogan, erdogan.eren@gmail.com
        %------------------------------------------------------------
        
        % Time Conversions
          % compute terrastrial time TT
            TT=h+(TAI_UTC+32.184)/3600;
          % compute terrastrial time UT1
            UT1=h-(UT1_UTC)/3600;            
        % Time since epoch J2000 in julian centruies
          [mjd,jd,t] = this.julianDay(y,m,d,TT);
        % Parameters of precession angel
          [za,thetaa,zetaa] =this. precession(t);
        % Parameters of nutation angel 
          [eps0,deps,dpsi] = this.nutation(t);  
        % Parameters of earth rotation angel 
          [gast] = this.jul2gast(UT1,t);
        % Transformation matrixes        
          % Precession,P
            P = this.rot('z',-za)*this.rot('y',thetaa)*this.rot('z',-zetaa);
          % Nutation, N
            N=this.rot('x',-eps0-deps)*this.rot('z',-dpsi)*this.rot('x',eps0);
          % Earth rotation, Er
            gast_deg=gast*360/24; % conversion from hour to degree
            Er=this.rot('z',gast_deg);
          % Polar motion, Pm
            sec2rad=1/3600*(pi/180); % second to radian
            xp=xp*sec2rad;
            yp=yp*sec2rad;
            Pm=[1   0   xp
                0   1   -yp
                -xp yp  1];
        % Rototation matrix from ICRS to ITRF
          E=Pm*Er*N*P;
        % First order time derivative of rotation matrix
          S=[0  1 0
             -1 0 0
             0  0 0];
          dEr_dt=w_ecef*S*Er;
          dE=Pm*dEr_dt*N*P;
          
          
        end
        function [gast,gmst] = jul2gast(this,ut1,t)
        % JUL2GAST returns the Greenwich Actual (true) Siderial Time (GAST)
        %
        % HOW gast = jul2gast(ut1,t)
        % IN  t    - time since epoch J2000.0        [in Julian centuries 
        %                                             (of 36525 days)]
        %     ut1  - universal time Greenwich        [decimal hours]
        % OUT gast - Greenwich Actual Siderial Time  [decimal hours !!!]
        %
        % NB  Due to the accuracy level of NUTATION, the result will be 
        % accurate to about 0.1 s (of time).

        %------------------------------------------------------------------
        % Sneeuw/Zebhauser, IAPG, TU Munich                    04/01/96
        %------------------------------------------------------------------
        % uses NUTATION
        %
        %------------------------------------------------------------------
        % revision history
        % .Oct.2001,NS: translation of original PREZWINK.M 
        %------------------------------------------------------------------
        % remarks
        %
        %------------------------------------------------------------------

        % Here we go
        gmst = (ut1*3600 + 24110.54841 + 8640184.812866*t + 0.093104*t.^2 ...
                - 6.2e-6*t.^3) /3600;
        [eps0,deleps,delpsi] = this.nutation(t);
        epsi = (eps0 + deleps) / 180 * pi;
        eq   = delpsi .* cos(epsi);
        gast = gmst + eq/15;
        gast = rem(rem(gast,24)+24,24);
        gmst = rem(rem(gmst,24)+24,24);

        end     
        function [mjd,jd,t] = julianDay(this,y,m,d,h)

        % HOW [t,mjd] = julian2000(y,m,d,ut1)
        % IN  y   - year (4-digit!)
        %     m   - month
        %     d   - day
        %     ut1 - universal time Greenwich      [decimal hours]
        % OUT t   - time since the epoch J2000.0  [in Julian centuries (of 36525 days)]
        %     mjd - modified Julian day           [days]
        %
        % NB  The transformation is exact.

        %-----------------------------------------------------------
        % Sneeuw/Zebhauser, IAPG, TU Munich                  04/01/96
        %-----------------------------------------------------------
        % uses none
        %
        %------------------------------------------------------------------
        % revision history
        % .980121NS: y/m/d-tests adapted. Vectorial input of y/m/d allowed.
        % .980122NS: MJD-output.
        % .Oct.2001,NS: translation of original NUTWINK.M 
        %------------------------------------------------------------------
        % remarks
        %
        %------------------------------------------------------------------
        % modifications
        % ut1-->h
        % function name : JULIAN2000 -->julianDay
        % error checking
        if any(m(:)>12 | m(:)<1), error('Month between 1 and 12'), end
        if any(d(:)>31 | d(:)<1), error('Day between 1 and 31'), end
        if any(rem(y(:),1) ~=0) | any(rem(m(:),1) ~=0) | any(rem(d(:),1) ~=0) 
           error('Year, Month, Day must be integer')
        end

        % The auxiliary variable is the Julian Day
        jd  = 367*y - floor(7*(y+floor((m+9)/12))/4);
        jd  = jd + floor(275*m/9) + d + 1721014 + h/24 - 0.5;
        t   = (jd-2451545)/36525;
        mjd =  jd-2400000.5;				% modified jd
        end   
        function [year, month, day, hour, minute, second] = jd2date(this,jd)
        %JD2DATE Gregorian calendar date from modified Julian day number.
        %
        %   [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND] = JD2DATE(JD) returns the
        %   Gregorian calendar date (year, month, day, hour, minute, and second)
        %   corresponding to the Julian day number JD.
        %
        %   Start of the JD (Julian day) count is from 0 at 12 noon 1 JAN -4712
        %   (4713 BC), Julian proleptic calendar.  Note that this day count conforms
        %   with the astronomical convention starting the day at noon, in contrast
        %   with the civil practice where the day starts with midnight.
        %
        %   Astronomers have used the Julian period to assign a unique number to
        %   every day since 1 January 4713 BC.  This is the so-called Julian Day
        %   (JD). JD 0 designates the 24 hours from noon UTC on 1 January 4713 BC
        %   (Julian calendar) to noon UTC on 2 January 4713 BC.

        %   Sources:  - http://tycho.usno.navy.mil/mjd.html
        %             - The Calendar FAQ (http://www.faqs.org)

        %   Author:      Peter John Acklam
        %   Time-stamp:  2002-05-24 15:24:45 +0200
        %   E-mail:      pjacklam@online.no
        %   URL:         http://home.online.no/~pjacklam


           % Adding 0.5 to JD and taking FLOOR ensures that the date is correct.
           % Here are some sample values:
           %
           %  MJD     Date       Time
           %  -1.00 = 1858-11-16 00:00 (not 1858-11-15 24:00!)
           %  -0.75 = 1858-11-16 06:00
           %  -0.50 = 1858-11-16 12:00
           %  -0.25 = 1858-11-16 18:00
           %   0.00 = 1858-11-17 00:00 (not 1858-11-16 24:00!)
           %  +0.25 = 1858-11-17 06:00
           %  +0.50 = 1858-11-17 12:00
           %  +0.75 = 1858-11-17 18:00
           %  +1.00 = 1858-11-18 00:00 (not 1858-11-17 24:00!)

           ijd = floor(jd + 0.5);               % integer part

           if nargout > 3
              fjd = jd - ijd + 0.5;             % fraction part
              [hour, minute, second] = this.days2hms(fjd);
           end

           % The following algorithm is from the Calendar FAQ.

           a = ijd + 32044;
           b = floor((4 * a + 3) / 146097);
           c = a - floor((b * 146097) / 4);

           d = floor((4 * c + 3) / 1461);
           e = c - floor((1461 * d) / 4);
           m = floor((5 * e + 2) / 153);

           day   = e - floor((153 * m + 2) / 5) + 1;
           month = m + 3 - 12 * floor(m / 10);
           year  = b * 100 + d - 4800 + floor(m / 10);
        end     
        function [hour, minute, second] = days2hms(this,days)
        %DAYS2HMS Convert days into hours, minutes, and seconds.
        %
        %   [HOUR, MINUTE, SECOND] = DAYS2HMS(DAYS) converts the number of days to
        %   hours, minutes, and seconds.
        %
        %   The following holds (to within rounding precision):
        %
        %     DAYS = HOUR / 24 + MINUTE / (24 * 60) + SECOND / (24 * 60 * 60)
        %          = (HOUR + (MINUTE + SECOND / 60) / 60) / 24

        %   Author:      Peter John Acklam
        %   Time-stamp:  2002-03-03 12:52:02 +0100
        %   E-mail:      pjacklam@online.no
        %   URL:         http://home.online.no/~pjacklam


           second = 86400 * days;
           hour   = fix(second/3600);           % get number of hours
           second = second - 3600*hour;         % remove the hours
           minute = fix(second/60);             % get number of minutes
           second = second - 60*minute;         % remove the minutes
        end        
        function [R]=rot(this,axis,angel)
        % FUNCTION
        %   Euler in counterclockwise direction around given axis
        % INPUTS
        %   angel : rotation angel in degree
        %   axis  : rotation axis. Takes 'x','y' or 'z'
        % OUTPUTS
        %   R     : rotation matrix
        %--------------------------------------------------------
        if strcmp(axis,'x')
            R=[1  0            0
               0  cosd(angel)  sind(angel)
               0  -sind(angel) cosd(angel)];
        elseif strcmp(axis,'y')
            R=[cosd(angel)   0    -sind(angel)
               0             1    0
               sind(angel)   0    cosd(angel)];
            
        elseif strcmp(axis,'z')
            R=[cosd(angel)   sind(angel)    0
               -sind(angel)  cosd(angel)    0
               0             0              1];            
            
        else
            error('wrong axis identifier')
        end
        end
        function [R]=orb2geo(this,x,y,z,vx,vy,vz)
        %*******************************************************
        %FUNCTION
        %   Transformation matrix from orbital coordinate system given by
        %   radial(eR),along-track(eT) and cross(eN)directions to geocentric
        %   reference system
        %INPUTS
        %   x,y,z   : position component of satellite in geocentric system
        %   vx,vy,vz   : velocity component of satellite in geocentric system
        %OUTPUTS
        %   ax ay az: acceleration global cartesian system (XYZ)
        %REFERANCE
        %   Montenbruck and Gill, Satellite Orbits
        %--------------------------------------------------------------
        
        r=[x;y;z]; v=[vx;vy;vz];
        eR=-r/(norm(r));  eN=-cross(r,v)/norm(cross(r,v)); eT=cross(eN,eR);
        R=[eR eT eN];        
        
        end
   

    end
    
end

