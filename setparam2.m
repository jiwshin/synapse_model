function setparam2(pw)

global Aca Alocal cb  kcloc kc1 kc2 Clow 
global Nsite Nmax Smax bSyt fw bw syt0
global ksb Ksmx Kds hs ksr pnb
global knb Knmx Kdn hn knr
global kinit kz kzr ke z1eff sigma kbeta
global g h n0 s0 x0 betamin betaSlope gamma Vstep;

    cb = 0.05;    

    idx = 1;
    Aca  = pw(idx); idx=idx+1;
    Alocal  = pw(idx); idx=idx+1;
    kcloc = pw(idx)*1e-3; idx=idx+1; 
    kc1 = pw(idx)*1e-3; idx=idx+1; %1/ms
    kc2 = pw(idx)*1e-3; idx=idx+1; %1/ms
    Clow = pw(idx); idx=idx+1;
    Nsite = pw(idx); idx=idx+1; %number of AZ
    Nmax = pw(idx); idx=idx+1; %number of SV docking site per AZ

    x0  = pw(idx); idx=idx+1; 
      
    idx=idx+1;
    knb = pw(idx)*1e-3; idx=idx+1; %1/ms
    knr = pw(idx)*1e-3; idx=idx+1;
    Knmx = pw(idx)*1e-3; idx=idx+1; %1/ms
    if isnan(Knmx)
        Kdn = pw(idx)*1e-3; idx=idx+1;
    else
        Kdn = pw(idx); idx=idx+1;
    end
    hn = pw(idx); idx=idx+1; 
    bSyt = pw(idx); idx=idx+1; 
    Vstep = pw(idx); idx=idx+1; 

    idx=idx+1;
    ke  = pw(idx)*1e-3; idx=idx+1;
    z1eff = pw(idx); idx=idx+1; 

    idx=idx+1;
    betamin= pw(idx)*1e-3; idx=idx+1;  
    kbeta = pw(idx)*1e-3; idx=idx+1; 
    pnb = pw(idx); idx=idx+1;
    
    n0 = Nmax*knb/(knb+knr);
    g = n0/Nmax;
    
    kinit = [cb, cb, knb, knb];    %[cb, clocal, knb, knglobal]
    fprintf('n0 = %1.1f \n', n0);


%%
    if bSyt
        Syts = Knmx-knb;
        %Syt7 (Table 3, Brandt, ACS, 2012)
        kon = 7*1e-3; %1/uM/ms
        koff = 10*1e-3; %1/ms
            %For Syt1
            %kon = 15*1e-3; %1/uM/ms
            %koff = 600*1e-3; %1/ms
        bs = bSyt;
        fw = [Vstep-1:-1:1].*kon;
        bw = ones(1, Vstep-1);
        for idx=1:Vstep-1
            bw(idx) = idx*bs^(idx-1);
        end
        bw = bw.*koff;
        %bw = [1, 2*bs, 3*bs^2, 4*bs^3].*koff; 
        syt0 = initSyt(cb, fw, bw).*Syts;
    else
        Vstep = 0;
    end

%************************************************************
    %{
    n0 = Nmax*knb / (knr + knb*(1 + ksb/ksr)) ;
    h = n0/Nmax;
    s0 = ksb*n0/(ksr + ksb*n0/m0);
    g = s0/Nmax;
    %}
    
    %{
    g0 = ksb/ksr * Smax/Nmax
    a0 = (g0-1)*(knb+knr)
    b0 = Nmax*((knb+knr)-(g0-1)*knb)
    c0 = -knb*Nmax*Nmax
    if g0 == 1
        n0 = knb*Nmax / (knb+knr)
    else
        n0 = roots([a0 b0 c0])
        n0 = min(abs(n0))
    end
    s0 = g0*n0*Smax / (g0*n0+Nmax-n0)
    %}
    
    % set en=1 (only in reverse)
    %{
    a0 = ksb*(knb+knr);
    b0 = Nmax/Smax*ksr*(knb+knr) - Nmax*knb*ksb;
    c0 = -knb*ksr*Nmax*Nmax/Smax;
    
    n0 = roots([a0 b0 c0]);
    n0 = max(n0);
        
    s0 = n0*ksb*Smax/(Nmax/Smax*ksr + n0*ksb);
    %}
    
    %{
    % rdx{ves} =  rxn involving ves 
    rdxv = cell(5,1);
    rdxv{1} = [1, 4, 7, 10];
    rdxv{2} = [1, 2, 7, 8];
    rdxv{3} = [2, 8];
    rdxv{4} = [4, 10];
    rdxv{5} = [4, 5, 10];

    % Qv = zeros(5, 1);
    % for idx=1:5
    %     Qv(idx) =  sum(rate(rdxv{idx}+1));
    % end

    % vdxr{rxn} = ves affected by rxn
    vdxr = cell(10,1);
    vdxr{1} = [1, 2];
    vdxr{2} = [2, 3];
    vdxr{4} = [1, 4, 5];
    vdxr{5} = [5];
    vdxr{7} = [1,2];
    vdxr{8} = [2,3];
    vdxr{10} = [1, 4, 5];
    %}

        
    %{
    Rn = kncw/knr;
    Rs = kscw/ksr;
    
    denom = 1 + Rn + Rn.*Rs;
    ninfw = Nmax*Rn./denom;
    sinfw =  Nmax*Rn.*Rs ./ denom;
    
    figure(23); plot(cw, ninfw, cw, sinfw); title('ninf sinf');

    %}
        
    %fprintf('ksb = %1.2f\n', kscw(1)*1000);

end