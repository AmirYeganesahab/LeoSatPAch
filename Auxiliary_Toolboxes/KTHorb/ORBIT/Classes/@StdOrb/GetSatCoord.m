function ret = GetSatCoord(std, prn, t)

%GetSatCoord - computes satellite's coordinates and clock correction for a
%given epoch, using standard ephemeris
%ret = GetSatCoord(std, prn, t)
%std - standard orbit in polynom form (class StdOrb) created by CompStdOrb
%prn - satellite number
%t - eopch in [FateTime]

%Written by Milan Horemuz, last modified 2004-11-01


if (t - std.dTbeg) < 0  | (t - std.dTend) > 0
    t
    error('Given epoch is outside of the standard ephemeris');
end
en = FindEn(std.dTimeSet, t);  %Finds epoch number en
if en < 1
    t
    error('Should not be here: something is wrong in findEn');
end
ent = FindEn(std.dtClTime, t);
if ent < 1
    t
    error('Should not be here: something is wrong in findEn finding dtClTime');
end

ret = struct('x',0,'y',0,'z',0,'dt',0, 'Tgd',0);
ret.Tgd = std.dTGD(prn, en);
coefX = zeros(1,size(std.dXCoef,3));
coefY = coefX;
coefZ = coefX;
coefDT =coefX;
Xmu = zeros(1,2);
Ymu = Xmu;
Zmu = Xmu;
%Tmu = Xmu;
coefX(:) = std.dXCoef(prn, en, :);
if abs(coefX(length(coefX))) < 1e-12
    prn
    error('No data for given prn');
end
coefY(:) = std.dYCoef(prn, en, :);
coefZ(:) = std.dZCoef(prn, en, :);
%coefDT(:) = std.dDtsCoef(prn, en, :);
Xmu(:) = std.dXmu(prn, en, :);
Ymu(:) = std.dYmu(prn, en, :);
Zmu(:) = std.dZmu(prn, en, :);
%Tmu(:) = std.dTmu(prn, en, :);
tt = t - std.dTimeSet(en);
ret.x = polyval(coefX, tt, [], Xmu);
ret.y = polyval(coefY, tt, [], Ymu);
ret.z = polyval(coefZ, tt, [], Zmu);
%ret.dt = polyval(coefDT, tt, [], Tmu);
%coefVX = polyder(coefX);
%tv = (tt - Xmu(1))/Xmu(2);
%ret.vx = polyval(coefVX, tv);

coefVX = coefX;
coefVY = coefY;
coefVZ = coefZ;
n = length(coefX);
coefVX(n) = []; %delete the last coefficient (first derivative = 0)
coefVY(n) = [];
coefVZ(n) = [];
n = n-1; %length of coefVX
for i=1:n
    coefVX(i) = coefVX(i)*(n-i+1);
    coefVY(i) = coefVY(i)*(n-i+1);
    coefVZ(i) = coefVZ(i)*(n-i+1);
end
tv = (tt - Xmu(1))/Xmu(2);
ret.vx = polyval(coefVX, tv)/Xmu(2);
tv = (tt - Ymu(1))/Ymu(2);
ret.vy = polyval(coefVY, tv)/Ymu(2);
tv = (tt - Zmu(1))/Zmu(2);
ret.vz = polyval(coefVZ, tv)/Zmu(2);

targ = std.dtClTime(ent:ent+1) - std.dtClTime(ent);
[coefDT, S, mu] = polyfit(targ, std.dClErr(prn, ent:ent+1), 1);
tt = t - std.dtClTime(ent);
ret.dt = polyval(coefDT, tt, [], mu);
