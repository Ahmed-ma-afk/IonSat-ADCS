p = fipref;
p.NumericTypeDisplay = 'short';
p.FimathDisplay = 'none';
p.LoggingMode = 'on';
F = fimath('OverflowAction','Wrap',...
     'RoundingMethod','Floor',...
     'CastBeforeSum',false);
warning off
format compact

a = fi(pi, true, 8, 0.1, 0)

b = fi(exp(1), true, 8, 0.1, 0)

F.ProductMode = 'FullPrecision';
F.SumMode = 'FullPrecision';
a.fimath = F;
b.fimath = F;

a

b

a*b

a+b