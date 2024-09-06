function IonoPar = getIonoPar(brObs)
% This function returns the broadcast Ionospheric parameters from a 
% navigation file header

IonoPar = [brObs.dAlfaIon1,brObs.dAlfaIon2,brObs.dAlfaIon3,brObs.dAlfaIon4;
          brObs.dBetaIon1,brObs.dBetaIon2,brObs.dBetaIon3,brObs.dBetaIon4];