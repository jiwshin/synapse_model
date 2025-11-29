%inherited from ni7, but gaussian local Ca inc is implemeted

function [vout, vw, kw] = ni7_2(vin, isi)

global dt Nsite T bAP
%global ksb Ksmx Kds hs ksr
global knb Knmx  knr gamma kzr
global kz ke z1eff n0 x0 Smax
global Aca  Alocal cb  kc1 kc2 Clow
global kinit kstat Nmax
global  bSyt Vstep syt0 sytstat

    if isempty(kstat)
        kstat = kinit;  %[cb, cloc, knb, ksb, knb, ksb];
    end

    if bSyt 
        fwloc =  0.2; %0.189; %full width
        tloc = 0.25; %timing of peak local Ca
        sdloc = fwloc/2.355; 
        Aloc = Alocal/normpdf(tloc, tloc, sdloc); % normpdf(t, m, sd)
        if isempty(sytstat)
            sytstat = syt0; 
        end
    end
    
    % vin = vdx by site matrix {e n}
    v = vin;
    
    nisi = ceil(isi/dt);
    kw = zeros(nisi, 3);
    idx = 1;    
    kw(idx,:) = kstat(1:3); %[c, cloc, kn]    
    vw = zeros(nisi, 2);
    vw(idx,:) = sum(vin,2)';
    if bSyt
        sytw = zeros(nisi, Vstep);
        sytw(idx,:) = sytstat;        
    end
    
    %T = T +dt;  
 
    c_ = kstat(1) - cb;
    %cloc = kstat(2) - cb;
    knglobal = kstat(4);

    t = 0;

    while (t<isi)
        if bAP && t==0
            c_ = c_  + Aca;
            dt_ = 0.02;
        else
            dt_ = dt;
        end

        if bSyt>0
                t_ = t ;
                while(t_ <=  t+dt)
                    cloc = c_ + Aloc*normpdf(t_, tloc, sdloc);
                    syt = sytw(idx,:);
                    syt = syt + deriv(syt, cloc, dt_);
                    t_ = t_ + dt_;
                end    
                knglobal = knb + Knmx*syt(Vstep);
                knw = knglobal*ones(1,Nsite);               
        end        

        %v(5) = z1, and z1 lowers kn or ks
         switch z1eff
            case 1
                knw = knw.*exp(-v(3,:));                
            case 2 
                knw = knw.*exp(-1*v(3, :)/n0);
            case 3
                knw = knw.*(1-v(3,:)/Nmax); 
        end
        
        for site = 1:Nsite
            kn = knw(site);

            % vin = vdx by site matrix {e  n }            
            transM = [-kn,      knr   ;...                        
                                kn,     -knr  ];

            v(:, site) = v(:, site) + transM*v(:, site)*dt;      
        end %for site          
        
        c1 = max(c_-Clow, 0);
        c2 = min(Clow, c_);
        dc = -(kc1*c1 + kc2*c2)*dt;
        c_ = c_ + dc;

        t = t+dt;    
        T = T + dt;
        idx = idx+1;

        kw(idx,:) = [c_+cb, cloc + cb, mean(knw)];        
        if bSyt
            sytw(idx,:) = syt;
        end
        vw(idx,:) = sum(v, 2)';
    end
    vout = v;
    kstat = [kw(end,:), knglobal];            
    if bSyt
        sytstat = sytw(end,:);
    end

end

function [dvdt] = deriv(syt, ca, dt)
	global fw  bw Vstep
	    fwc = fw.*ca;
        forw = zeros(1, Vstep);
        back = zeros(1, Vstep);
        dvdt = zeros(1, Vstep);

        forw(1, 1:Vstep-1) = fwc.*syt(1, 1:Vstep-1);
        for idx=2:Vstep
            back(idx) = bw(idx-1)*syt(idx);
        end
		
		dvdt(1) = back(2) - forw(1);

        for idx=2:Vstep-1
            dvdt(idx) = forw(idx-1) + back(idx+1) - back(idx) - forw(idx);
        end

        dvdt(Vstep) = forw(Vstep-1) - back(Vstep);
		dvdt = dvdt .*dt;
end
