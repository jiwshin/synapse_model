function [initv] = initSyt(c, fw, bw)
    global Vstep
    fwc = fw.*c;
	fbw = fwc./bw;
    initv = zeros(1, Vstep);
% Syts = v1+v2+v3+v4+v5
    denom = 1;
    add = 1;
    for idx=1:Vstep-1
        add = add*fbw(idx);
        denom = denom + add;
    end
    %	initv(1) = 1 / (1 + fbw(1) + fbw(1)*fbw(2) + fbw(1)*fbw(2)*fbw(3) + fbw(1)*fbw(2)*fbw(3)*fbw(4));
    initv(1) = 1/denom;
    for idx=2:Vstep
        initv(idx) = fbw(idx-1)*initv(idx-1);
    end
%  	initv(2) = fbw(1)* initv(1); 
%  	initv(3) = fbw(2) * initv(2);
%  	initv(4) = fbw(3) * initv(3);
%  	initv(5) = fbw(4) *initv(4); 

end