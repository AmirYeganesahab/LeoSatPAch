function vpeEph = RemoveEpoch(vpeEphi,tRef)

%function vpeEph = RemoveEpoch(vpeEphi,tRef)
%Removes epoch tRef (FateTime) from vpeEphi (CoordEph)
%Returns vpeEph (CoordEph) with removed epoch

%Written by Milan Horemuz, last modified 2005-02-04


vpeEph = vpeEphi;
t = get(vpeEph,'dTime');  %extract vector of eopchs
nep = length(t);  %number of epochs in vpeEph
epnr = -1;
for i = 1:nep
    if abs(t(i) - tRef) < 1e-5
        epnr = i;
        break;
    end
end
if epnr < 0
    tRef
    error('Could not find this epoch in vpeEph');
end

vpeEph.dX(:,epnr) = [];
vpeEph.dY(:,epnr) = [];
vpeEph.dZ(:,epnr) = [];
vpeEph.dDts(:,epnr) = [];
vpeEph.iPRN(:,epnr) = [];
vpeEph.dTime(:,epnr) = [];
