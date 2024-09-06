function A = ReadEphPrecise(varargin)
%function A = ReadEphPrecise(A,sFileName_1 ... sFileName_n)

%ReadEphPrecise reads precise ephemeris (coordinates and clock corrections) from sp3
%file
%A = ReadEphPrecise(A,sFileName) 
%A instance of class CoordEph, which is container of coordinate ephemeris
%sFileName 1 ... n - file(s) containing sp3 ephemeris
%A is to be used in CompStdOrb

%Written by Milan Horemuz, last modified 2005-02-03

A = CoordEph; %create an instance
i = 1; %counter of epochs for which the satellite coordinates are listed in PE
for npe = 2:nargin
	%read header
	first ='1'; % variable for testing end of  header
	fid=fopen(varargin{npe}); %open file
	if fid < 1
        fprintf('Could not open file %s', sFileName);
        return;
	end
	tline = fgetl(fid);
	tline(1:3) = [];  %delete the first 3 characters 
	epoch = sscanf(tline,'%f');
	Time0 = FateTime(epoch);
	tline = fgetl(fid); %skip 2 tline
	tline = fgetl(fid);
	tline(1) = [];  %delete + character
	MofSat = sscanf(tline, '%i'); %read number of satellites
	while first ~= '*'
        tline = fgetl(fid);
        if length(tline) < 3
            fprintf('Not a valid sp3 file');
            return;
        end
        first = tline(1);
	end
	first ='1'; % variable for testing the line with time
	while feof(fid) == 0 & tline(1) ~= 'E';
        tline(1) = [];  %delete * character
        epoch = sscanf(tline,'%f');
        A.dTime(i) = FateTime(epoch(1), epoch(2), epoch(3), epoch(4), epoch(5), epoch(6));  %date2jd(epoch(1), epoch(2), epoch(3), epoch(4), epoch(5), epoch(6));
        tline = fgetl(fid);
        while first ~= '*' 
            if tline(2) == 'R' %read in only GPS satellites
                tline = fgetl(fid);
                first = tline(1);
                if first == 'E'
                    break;
                end
                continue;
            end
            tline(1:2) = []; %delete PG
            prn = sscanf(tline(1:2), '%d'); %satellite number
            epoch = sscanf(tline,'%f'); 
            A.iPRN(prn,i) = epoch(1);  
            A.dX(prn,i) = epoch(2)*1000; %convert to [m]
            A.dY(prn,i) = epoch(3)*1000; 
            A.dZ(prn,i) = epoch(4)*1000; 
            A.dDts(prn,i) = epoch(5)/1e6; %convert to [s]
            tline = fgetl(fid);
            if tline(1) == 'E'
                break;
            end
            first = tline(1);
        end
        i = i+1;
        first ='1';
	end
	fclose(fid);
end
