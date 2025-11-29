%FWHM = 2.355 sigma = 0.189
isi = 200;
dt = 1;
dt_ = 0.02; 
T = 0;

Aca = 1;
kc1 =20e-3;
kc2 = 20e-3;
Clow = 0;

fwloc =  0.2; %0.189; %full width
tloc = 0.25; %timing of peak local Ca
pkcloc = 40; % peak loca Ca
sdloc = fwloc/2.355; 
Aloc = pkcloc/normpdf(tloc, tloc, sdloc); % normpdf(t, m, sd)

Npulse = 20;
nt = isi*Npulse;
ntloc = (isi+dt/dt_)*Npulse;
calocw = zeros(ntloc, 1); 
caw = zeros(nt, 1); 
tw = zeros(ntloc, 1); 
t=0;
c_ = 0;
idx = 1;
jdx = 1;
for pulse = 1:Npulse
    t=0;
     while (t<isi)
            if t==0
                c_ = c_  + Aca;
                dt_ = 0.02;
            else
                dt_ = dt;
            end
    
            t_ = t ;
            while(t_ <=  t+dt)
                cloc = c_ + Aloc*normpdf(t_, tloc, sdloc);
                t_ = t_ + dt_;
                 idx = idx+1;
                 tw(idx) = isi*(pulse-1)+t_;
                calocw(idx) = cloc;
            end    
            c1 = max(c_-Clow, 0);
            c2 = min(Clow, c_);
            dc = -(kc1*c1 + kc2*c2)*dt;
            c_ = c_ + dc;
            caw(jdx) = c_;
            jdx = jdx + 1;
    
            t = t +dt;    
            T = T + dt;
     end
end

% tw= 0:dt_:2; 
% caw = Aloc*normpdf(tw, tloc, sdloc);
figure(12); clf; plot(tw, calocw); hold on; plot(caw);