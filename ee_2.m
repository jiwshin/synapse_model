%inherited from ee.m but it uses ni7_2
% If bSyt =1, this code simulate STP assuming that kn is proportional to Syt7(Vstep)
% bsyt = cooperativity of calcium binding
% Vstep = how many Ca ions are bound to Syt7
% Cascade Ca binding to syt7 is implemented in deriv() of ni7.m
% [Syt7] = Knmx - knb
% kon and koff are defined in setparam.m


set(0, 'DefaultFigureWindowStyle', 'docked');
%clrcode = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];

global cb isi dt T bAP bSyt syt0 sytstat
global pnb Nsite Nmax g n0
global vinit kinit kstat 
%{
global cw kncw kscw betaw gamw
global knb Knmx Kdn hn ksb Ksmx Kds hs
global betamin betaSlope gamma
%}

pfname = "param2.xls";
bAPw = [1 0];

if ~exist('rcvdata', 'var')    
    rcvdata = readtable('rcvwtee.xlsx');
end

if ~exist('stpdata', 'var')
    stpdata = readtable('stpwtee.csv');
end

if ~exist('bPlot', 'var')
    bPlot = 1;
end

if ~exist('pset', 'var')
    pset = readtable(pfname);
end
 
pw=pset.pw3; %pw1,2,3 (WT), pw4 (KD)
setparam2(pw);  % here kinit is set to [cb, cloc, knb, knglobal]

%%
%writetable(pset, pfname);

%% params 

dt = 1;
%isiw = 200;  clrcode = "#0072BD";
%isiw = 25;  clrcode = "#7E2F8E";
%isiw = [200 25];  clrcode = ["#0072BD",  "#7E2F8E"];
isiw = [200, 100, 50, 25];  clrcode = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E"];
Nf = length(isiw);   %num of frq

Ntrial = 10;

Thead = -10;
Ttail = 10;
Npuls = 20; 
cb = 0.05;

rls = cell(Nf, 1); 
nz = zeros(Npuls, Nf);
z  = zeros(Npuls, Nf);
vrec = cell(Ntrial, Nf);
krec = cell(Ntrial, Nf);
v_kin = cell(Nf, 1);
k_kin = cell(Nf, 1);
c = cb;

% vinit = vdx by Nsite matrix
% kinit = [c, clocal, kn, ks];

%% Main
tic
for fdx = 1:Nf
    isi = isiw(fdx);
    T = 0;
    c = cb;
    rls{fdx} = zeros(Npuls, Ntrial);
    
    for ndx=1:Ntrial

        % initialize v -------------------------------------- row of vinit is [en, es, n,s, z1, z2, x]
        vinit = zeros(2, Nsite);

        for site = 1:Nsite
            for nn = 1:Nmax
                 if rand() < g     % g = ninf/Nmax
                    vinit(2, site) = vinit(2, site) + 1; %n
                 end
            end
        end
        vinit(1, :) = Nmax - vinit(2, :); % e
        % ____________________________________
        
        %Header
        vrec{ndx, fdx} = [vrec{ndx, fdx}; repmat(sum(vinit, 2)', abs(Thead), 1)];
        krec{ndx, fdx} = [krec{ndx, fdx}; repmat(kinit(1:3), abs(Thead), 1)];
        if bSyt
            
        end

        %Train --------------------------------------------
        bAP = bAPw(1); %boolean AP
        v = vinit;
        kstat = kinit;
        sytstat = syt0;

        for idx = 1:Npuls            
            for site=1:Nsite
                n = v(2, site);

                if bAP
                    zn = pnb*n;
                else
                    zn = 0;
                end

                v(1,site) = v(1,site) + zn;  % e
                v(2,site) = v(2,site) - zn;  % n
                    
                 rls{fdx}(idx, ndx) = rls{fdx}(idx, ndx) + zn;

            end %for site

            if bSyt
                [vout, vw, kw] = ni7_2(v, isi);  %Syt7
            else
                [vout, vw, kw] = ni(v, isi);  
            end
            vrec{ndx, fdx} = [vrec{ndx, fdx}; vw];
            krec{ndx, fdx} = [krec{ndx, fdx}; kw];        

            v = vout;
            %kstat = kw(end,:);
            T = T + dt;
            
        end %for idx = 1:Npuls __________________

        %Tail
            bAP = bAPw(2);
            if bSyt
                [vout, vw, kw] = ni7_2(v, Ttail);  % Syt7
            else
                [vout, vw, kw] = ni(v, Ttail);  % numerical integration
            end        
        
        vrec{ndx, fdx} = [vrec{ndx, fdx}; vw];
        krec{ndx, fdx} = [krec{ndx, fdx}; kw];

    end %for Ntrial
   
    v_kin{fdx} = cellmean_(vrec(:, fdx));
    k_kin{fdx} = cellmean_(krec(:, fdx));

     
     tx = Thead*dt:dt:T;    
     tx = tx(1:size(v_kin{fdx}, 1));
     
     if bPlot
         plotparams; 
     end    
    
    %rls{fdx} = Npuls by Ntrial
    z(:, fdx) = mean(rls{fdx}, 2);
    baseline = n0 * Nsite * pnb; 
    nz(:, fdx) = z(:, fdx)/baseline;
        
end %for Nf
toc

%% plot var

pulx = 1:Npuls;
pulx20 = 1:20;

rlsm_ = zeros(Npuls, Nf);
rlsvar_ = zeros(Npuls, Nf);
rlscumvar_ = zeros(Npuls, Nf);
rlscum_m_ = zeros(Npuls, Nf);
rlsm = zeros(Npuls, Nf);
rlsvar = zeros(Npuls, Nf);
rlscumvar = zeros(Npuls, Nf);
rlscum_m = zeros(Npuls, Nf);

if bSyt == 1
% Note rls = Npuls x Ntrial matrix    
    for fdx = 1:Nf
        [rlsvar_(:,fdx), rlsm_(:,fdx)] = var(rls{fdx}, 0, 2);
        rlscum = cumsum(rls{fdx}, 1);
        [rlscumvar_(:,fdx), rlscum_m_(:,fdx)] = var(rlscum, 0, 2);
    end

    figure(12); clf; 

    q = 1; %25; 
    q2 = q*q;

    for fdx = 1:Nf
        rlsm(:, fdx) = rlsm_(:, fdx).*q/Nsite;
        rlsvar(:, fdx) = rlsvar_(:, fdx).*q2/Nsite;
        plot(rlsm(:, fdx), rlsvar(:, fdx), 'o'); hold on
    end

    fitx=cat(1,rlsm(:, 1),rlsm(:, 2));
    fity=cat(1,rlsvar(:, 1),rlsvar(:, 2));

    userfit=fittype({'x^2','x'});
    f = fit(fitx, fity, userfit)

    xmax = floor(-1/f.a);
    x=linspace(0, xmax);
    y=(f.a)*x.^2 + (f.b)*x;
    plot(x,y,'k-'); hold on

%     rlsm_ = reshape(rlsm, Npuls*Nf, 1);
%     rlsvar_ = reshape(rlsvar, Npuls*Nf, 1);        
%     plot(rlsm_, rlsvar_, 'o'); hold on; 
    
    rlscum_m = rlscum_m_.*q/Nsite;
    rlscumvar = rlscumvar_.*q2/Nsite;
    plot(rlscum_m, rlscumvar, 'x')
    xlim([0 xmax*3]);  ylim([0 10]);
    datacursormode on; dcm = datacursormode(gcf); set(dcm,'DisplayStyle','window'); 
    dcm.UpdateFcn = @cursorinfo;    

    Nfit = 1/f.a;
    qfit = 1/f.b; 
    save('mc_', 'Nfit', 'qfit', 'rlsm', 'rlsvar', 'rlscum_m', 'rlscumvar', 'nz', 'Nsite');
end

%% Plot STP

figure(11); clf;
plot(pulx20, stpdata.f5, pulx20, stpdata.f10,  pulx20, stpdata.f20, pulx20, stpdata.f40);    
%plot(pulx20, stpdata.f5, pulx20, stpdata.f40);    
ylim([0 3]); hold on  
%title('nz, f5, f10, f20, f40'); 

if Nf == 1
    plot(pulx, nz(:,1), 'ok');     
else
    for idx = 1:Nf
        if mod(idx, 2) == 1
            plot(pulx, nz(:, idx), 'o', 'MarkerFaceColor', clrcode(idx)); hold on
        else
            plot(pulx, nz(:, idx), 'o',  'Color', clrcode(idx)); hold on        
        end
    end

end

%%
function txt = cursorinfo(~,info)
    idx = info.DataIndex;
    x = info.Position(1);
    y = info.Position(2);
    txt = [num2str(idx),  ',   ',  num2str(x) ',   ',  num2str(y)];    
end
