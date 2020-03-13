function U = potentialE(X,qdata,Iqdata,weight,sspFlag,XFix)
XAll = nan(1,length(sspFlag));
XAll(sspFlag) = X;
XAll(~sspFlag) = XFix;
Sigmavi = exp(XAll(end));
logp = -( log10(fcn_saxs_sphere_hmcmc(XAll(1:end-1),qdata))-log10(Iqdata) ).^2 / (2*Sigmavi) - log(sqrt(2*pi*Sigmavi));
U = -sum(weight.*logp);  
