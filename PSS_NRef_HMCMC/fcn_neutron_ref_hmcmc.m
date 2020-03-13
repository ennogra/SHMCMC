function varargout = fcn_neutron_ref_hmcmc(x1,qset)
global fitFlag x2 
global lbub
global lambda conv_method
global nsets
global interface_profile  nSlice nlayer
global binedp_delta_cutoff edp_0  edp
global plotEDPFlag

% --- Generate all the fitting parameters
x=zeros(length(fitFlag),1);
x(fitFlag==1) = x1;
x(fitFlag==0) = x2;
x = inversetransformx(x,lbub);

Inorm           = x(1);
Ibk             = x(2);
resdawave       = x(3);
resconstant     = x(4);
alphai_offset   = x(5);
sld_env         = x(6);
sld_sub         = x(7);
sigma_sub       = x(8);
layeredp        = reshape(x(9:end),nlayer,[]);

%% --- calculate and preview edp
rawedp = nan(2+nlayer,5);
rawedp(:,1:2) = [...
    0                       sld_env
    cumsum(layeredp(:,1))   layeredp(:,2)
    NaN                     sld_sub
    ];
rawedp(:,3) = 0;                                % add dispersion
rawedp(:,4) = [layeredp(:,3); sigma_sub; NaN];  % add roughness
rawedp(:,5) =  interface_profile;               % add interface profile
if interface_profile >= 2
    rawedp(:,6) = [layeredp(:,4);0;NaN];        % asymmetric factor
    rawedp(:,7) = rawedp(:,4);
end
edp_0 = edprofile(rawedp,nSlice);
edp = binedprofile(edp_0,binedp_delta_cutoff);

if plotEDPFlag == 1
    figure
    plot(edp_0(:,1),edp_0(:,2),'-');
    hold(gca,'on');
    plot(edp(:,1),edp(:,2),'o');
    %    assignin('base','edp_0',edp_0);
    xlabel(['depth (',char(197),')']);
    ylabel('SLD')
    z = edp(1,1);
    edp_step = rawedp(1,2);
    for ii=1:nlayer+1
        z = [z; rawedp(ii,1); rawedp(ii,1)];
        edp_step = [edp_step; rawedp(ii,2); rawedp(ii+1,2)];
    end
    z = [z; edp(end,1)];
    edp_step = [edp_step; rawedp(end,2)];
    plot(z,edp_step,'--b','linewidth',1.5);
    %     save('SLD.dat','edp_0','-ascii', '-tabs');
end

%% convert edp to delta in order to use x-ray reflectivigy code
edp(:,2) = edp(:,2)*lambda^2/2/pi;

%% calculate reflectivity
I_cal_set = cell(1,nsets);
qset_true = cell(1,nsets);
for iset = 1:nsets
    qz = qset{iset};
    alphai = alphai_offset + asind(qz*lambda/(4*pi));
    qz = 4*pi/lambda*sind(alphai);
    if resdawave == 0 && resconstant == 0 %--- no reoslution
        RR = parratt(qz,edp,lambda);
        I_cal = RR(:,2);
        I_cal = I_cal/Inorm + Ibk;
    else
        % --- with resolution
        resfactor = 3;
        nconvpoint = 11;
        RR_conv =  refconv(qz,edp,lambda,resdawave,resconstant,resfactor,nconvpoint,conv_method);
        I_cal = RR_conv(:,2);
        I_cal = I_cal/Inorm + Ibk;
    end
    I_cal_set{iset} = I_cal;
    qset_true{iset} = qz;
end

if nargout == 1
    varargout{1} = I_cal_set;
elseif nargout == 2
    varargout{1} = I_cal_set;
    varargout{2} = qset_true;
else
    error('Invalid number of output argument');
end