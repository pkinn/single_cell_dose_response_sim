function dr = returnDR(x,d)
%returnDR get dose response values from parameters and dose vector
%x should be [EC50, E0, Emax, HS]

dr = x(3)+ (x(2)-x(3))./(1+(d./x(1)).^x(4));

end

