function display(A)

fprintf('Year                     : %4d\n', A.year);
fprintf('Month                    : %4d\n', A.month);
fprintf('Day                      : %4d\n', A.day);
fprintf('Hour                     : %4d\n', A.hour);
fprintf('Minute                   : %4d\n', A.min);
fprintf('Seconds                  : %.5f\n', A.sec);
fprintf('GPS week                 : %4d\n', A.gweek);
fprintf('Seconds of GPS week      : %.5f\n', A.wsec);
fprintf('Day of week              : %4d\n', A.dweek);
fprintf('MJD                      : %.5f\n', A.MJD);
fprintf('Day of Year              : %.5f\n', A.DOY);
%struct('Year', A.year, 'Month', A.month, 'Day', A.day, 'Hour', A.hour, 'Minute', A.minute, 'Sec', A.sec, 'GPSweek', A.gweek, 'GPSsec', A.wsec);
