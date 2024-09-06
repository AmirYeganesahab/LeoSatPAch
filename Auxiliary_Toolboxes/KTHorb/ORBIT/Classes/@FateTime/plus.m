function r = plus(a,b)

if isa(a, 'FateTime') & isa(b, 'double')
    pom = a.wsec + b;
    r = FateTime(a.gweek, pom);
elseif isa(b, 'FateTime') & isa(a, 'double')
    pom = a + b.wsec;
    r = FateTime(b.gweek, pom);
else
    err = sprintf('Operator plus in FateTime does not allow arguments %s %s', class(a), class(b));
    error(err);
end