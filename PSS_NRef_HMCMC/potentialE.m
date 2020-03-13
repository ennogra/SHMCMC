function U = potentialE(X,nsets,qset,Iset,Iweight_set)
Sigmavi = exp(X(end-nsets+1:end));
I_cal_set = fcn_neutron_ref_hmcmc(X(1:end-nsets),qset);
U = 0;
for iset = 1:nsets
    I_cal = I_cal_set{iset};
    I = Iset{iset};
    Iweight = Iweight_set{iset};
    logp = -( log10(I_cal)-log10(I) ).^2 / (2*Sigmavi(iset)) - log(sqrt(2*pi*Sigmavi(iset)));
    U = U - sum(Iweight.*logp);
end