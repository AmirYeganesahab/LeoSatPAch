function A = MJD2date(A)

%MJD2date converts MJD to GPS time
%A = MJD2date(A)
%A is FateTime object

%Written by Milan Horemuz, last modified 2004-11-01

jd = A.MJD + 2400000.5;
week = (jd - 2444244.5)/7;
A.gweek = floor(week);
A.wsec = (week - A.gweek)*86400*7;
A = GPS2Date(A);



