    %global Nsite
    
    if fdx ==1
        figure(1); clf; title('e'); hold on
        figure(2); clf; title('n'); hold on
        figure(5); clf; title('c, clocal'); hold on
        figure(6); clf; title('kn'); hold on       

    end
    figure(1); plot(tx, v_kin{fdx}(:, 1)); hold on
    figure(2); plot(tx, v_kin{fdx}(:, 2)); hold on
    figure(5); plot(tx, k_kin{fdx}(:, 1)); hold on; plot(tx, k_kin{fdx}(:, 2));
    figure(6); plot(tx, k_kin{fdx}(:, 3)); hold on; 


  