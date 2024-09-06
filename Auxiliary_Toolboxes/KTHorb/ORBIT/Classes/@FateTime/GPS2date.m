function A = GPS2date(A)

%GPS2date converts datum to GPS time and MJD
%A = GPS2date(A)
%A - object of type FateTime

%Written by Milan Horemuz, last modified 2005-02-09


   jd = A.gweek*7 + A.wsec/86400 + 2444244.5;
   A.MJD = jd-2400000.5;
   a = floor(jd+0.5); 
   b = a + 1537; 
   c = floor((b-122.1)/365.25);
   d = floor(365.25*c); 
   e = floor((b-d)/30.6001);
   f = jd+0.5;
   A.day = b - d - floor(30.6001*e); % + (int) modf(jd+0.5,&pom);
   A.month = e-1-12* floor(e/14);
   A.year = c - 4715 - floor((7+ A.month)/10);
   A.dweek = floor(A.wsec/86400.);
   pom = A.wsec/3600 - A.dweek*24;
   A.hour = floor(pom);
   pom = (pom - A.hour)*60;
   %A.min = floor(pom);
%   A.sec = (pom - A.min)*60;
   A.sec = A.wsec - A.dweek*86400 - A.hour*3600;  % - A.min*60;
   A.min = floor(A.sec/60);
   A.sec = A.sec - A.min*60;
