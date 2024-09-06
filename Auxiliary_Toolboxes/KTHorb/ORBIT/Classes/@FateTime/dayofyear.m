function [doy, ep] = dayofyear(ep)

%computes day of year
%Written by Milan Horemuz 2005-03-04

%y = get(ep,'year');
ep0 = FateTime(ep.year,1,1,0,0,0);
%doy = floor(get(ep,'MJD') - get(ep0,'MJD'))+1;
ep.DOY = floor(ep.MJD - ep0.MJD) + 1;
doy = ep.DOY;
