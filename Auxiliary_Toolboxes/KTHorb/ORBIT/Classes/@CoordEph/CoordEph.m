function ret = CoordEph()

%CoordEph - Coordinate ephemeris container - cartesian coordinates and clock corrections
% Can contain coordinates from precise or broadcast ephemeris
%X-coordinate                   dX
%Y-coordinate                   dY
%Z-coordinate                   dZ
%Satellite clock correction     dDts
%Satellite nr                   iPRN
%Time                           dTime
% Group delay                   dTGD
%Written by Milan Horemuz, last modified 2005-01-31 by Johan Vium Andersson


jd = FateTime;
s = struct('dX', 0, 'dY', 0, 'dZ', 0, 'dDts', 0, 'iPRN',0, 'dTime', jd,'dTGD',0);
ret = class(s,'CoordEph');