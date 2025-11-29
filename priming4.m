global Syts dt kon koff fw bw Tmax Nstep Nmax

%% time prms
dt = 0.2; %/ms
Tpre = 10;
Npuls = 20;
isi = 25;
isin = floor(isi/dt); %ms
Tmax = Tpre+ isi*Npuls + 200; %ms

%% prms for Ca

bLocalCa = 1;
crest = 0.05; %uM
Aca = 1; %uM
cathr = 1;
kc1=  200e-3; %1/ms
kc2 = 20e-3; %1/ms

%% syt prms

Syts = 100;

Nstep = 5;
bSyt7 = 1;
if bSyt7 
    %Syt7 (Table 3, Brandt, ACS, 2012)
    kon = 7*1e-3; %1/uM/ms
    koff = 10*1e-3; %1/ms
else 
    %Syt1
    kon = 15*1e-3; %1/uM/ms
    koff = 600*1e-3; %1/ms
end

b = .4;

fw = [4, 3, 2, 1].*kon;
bw = [1, 2*b, 3*b^2, 4*b^3].*koff; 

%% initial setting
Nmax = 6;
Nt = (Tmax+Tpre)/dt+1;
xt = zeros(Nt, Nstep);
tw = zeros(Nt, 1);
cw = zeros(Nt, 1)+crest;
vinit = init(crest);	
xt(1,:) = vinit*Syts;
tw(1) = -Tpre;
ca = crest;
cathr = 1;
ti = -Tpre/dt;
idx  = 2;
ipul = 0;

while(ti<Nt)
    if mod(ti, isin) == 0 && ti > 0 && ipul < Npuls
        ca = ca + Aca;
        ipul = ipul + 1; 
    end
	tw(idx) = ti*dt;
    xt(idx, :) = xt(idx-1,:) + deriv(xt(idx-1, :), ca);
    if bLocalCa
        ca1 = max(0, ca - cathr);
        ca2 = min(ca, cathr);
        ca1 = ca1*(1-kc1*dt);
        ca2 = crest + (ca2-crest)*(1-kc2*dt);
        ca = ca1+ca2;
    else
        ca = crest + (ca-crest)*(1-kc2*dt);
    end
    cw(idx) = ca;
	ti = ti + 1;
    idx = idx+1;
end
%%
figure(1); clf; plot(vinit*Syts); hold on; plot(xt(end, :));
figure(2); clf; plot(tw, xt(:,Nstep), tw, xt(:,Nstep-2), tw, xt(:,1)); hold on 
yyaxis right; plot(tw, cw); xlim([-Tpre, Tmax]);
%%
function [initv] = init(c)
    global fw bw Nstep
    fwc = fw.*c;
	fbw = fwc./bw;
    initv = zeros(1, Nstep);
% RP = v1+v2+v3+v4+v5+v6
	initv(1) = 1 / (1 + fbw(1) + fbw(1)*fbw(2) + fbw(1)*fbw(2)*fbw(3) + fbw(1)*fbw(2)*fbw(3)*fbw(4));
 	initv(2) = fbw(1)* initv(1); 
 	initv(3) = fbw(2) * initv(2);
 	initv(4) = fbw(3) * initv(3);
 	initv(5) = fbw(4) *initv(4); 

end

function [dvdt] = deriv(ves, ca)
	global fw  bw dt Nstep Nmax
	    fwc = fw.*ca;
        forw = zeros(Nstep, 1);
        back = zeros(Nstep, 1);
		forw(1) = fwc(1)*ves(1);
		forw(2) = fwc(2)*ves(2);
		forw(3) = fwc(3)*ves(3);
		forw(4) = fwc(4)*ves(4);
		%forw(5) = fwc(5)*ves(5); %*(Nmax-ves(6))/Nmax;

		back(2) = bw(1)*ves(2);
		back(3) = bw(2)*ves(3);
		back(4) = bw(3)*ves(4);
		back(5) = bw(4)*ves(5);
		%back(6) = bw(5)*ves(6);
		
		dvdt(1) = back(2) - forw(1);
		dvdt(2) = forw(1) + back(3) - back(2) - forw(2);
		dvdt(3) = forw(2) + back(4) - back(3) - forw(3);
		dvdt(4) = forw(3) + back(5) - back(4) - forw(4);
        dvdt(5) = forw(4) - back(5);
% 		dvdt(5) = forw(4) + back(6) - back(5) - forw(5);
% 		dvdt(6) = forw(5) - back(6);
		dvdt = dvdt .*dt;
end
