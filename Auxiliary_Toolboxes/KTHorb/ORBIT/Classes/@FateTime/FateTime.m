function A = FateTime(varargin)  

%FateTime -  class that stores time and performes conversion between various format
%A = FateTime() - fills in zeros to all variables
%A - class FateTime
%A = FateTime(y, mo, d, ho, min, sec)
%A = FateTime([y mo d ho min sec])
%y, mo, d, ho, min, sec - date and time
%A = FateTime(gw, ws) 
%gs - GPS week
%ws - seconds of GPS week 
%A = FateTime(MJD) 
%MJD - modified julian date
%dweek - day of GPS week

%Written by Milan Horemuz, last modified 2004-11-01


na = length(varargin);
switch na
    case 0
        ret = struct('MJD', 0, 'gweek', 0, 'dweek', 0, 'wsec', 0, 'year', 0, 'month', 0, 'day', 0, 'hour', 0, 'min', 0, 'sec', 0,'DOY',0);
        A = class(ret, 'FateTime');
    case 6  %input in georgian datum and time of day
        [y, mo, d, ho, min, sec] = deal(varargin{:});
        ret = struct('MJD', 0, 'gweek', 0, 'dweek', 0, 'wsec', 0, 'year', y, 'month', mo, 'day', d, 'hour', ho, 'min', min, 'sec', sec,'DOY',0);
         A = class(ret, 'FateTime');
         A = date2GPS(A);
    case 2  %input in GPS week and wsec
        [gw, ws] = deal(varargin{:});
        ret = struct('MJD', 0, 'gweek', gw, 'dweek', 0,'wsec', ws, 'year', 0, 'month', 0, 'day', 0, 'hour', 0, 'min', 0, 'sec', 0,'DOY',0);
        A = class(ret, 'FateTime');
        A = GPS2date(A);
    case 1 
        n = length(varargin{1});
        if n == 1 %input in MJD
            mjd = varargin{1};
            ret = struct('MJD', mjd, 'gweek', 0, 'dweek', 0,'wsec', 0, 'year', 0, 'month', 0, 'day', 0, 'hour', 0, 'min', 0, 'sec', 0,'DOY',0);
            A = class(ret, 'FateTime');
            A = MJD2date(A);
        elseif n >= 6
            in=varargin{1};
            y = in(1);
            mo = in(2);
            d = in(3);
            ho = in(4);
            min = in(5);
            sec = in(6);
            ret = struct('MJD', 0, 'gweek', 0, 'dweek', 0, 'wsec', 0, 'year', y, 'month', mo, 'day', d, 'hour', ho, 'min', min, 'sec', sec,'DOY',0);
            A = class(ret, 'FateTime');
            A = date2GPS(A);
        end
        
     otherwise
         na
         error('FateTime: Incorrect number of arguments');
 end



