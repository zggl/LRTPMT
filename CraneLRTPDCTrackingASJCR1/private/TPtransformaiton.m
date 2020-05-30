function [S1,U1,maxerr,meanerr]=TPtransformaiton(gridsize1,lpvux,dep1,domain1,hull,svtol, keep,primebase,pbaseindx)
%% TP transformation, same as:
    %   [S U] = tptrans(lpv, dep1, domain, gridsize, 'close');
    % sampling
    p = primes(primebase);    
    p=p(pbaseindx:end); %  p_1=5 5x2
    TP.siz=prod(gridsize1);
    TP.lpv=lpvux;
    TP.dep=dep1;
    TP.domain =domain1;
    TP.gridsize=gridsize1;
    Datanum=TP.siz;
    Xij= DIST_HAMMERSLEYDiy(domain1, Datanum,p);
    TP.X_scaled=Xij;
    lpvdata1 = sampling_lpvud(TP);
    %lpvdata1 = sampling_lpv(lpvux, dep1, domain1, gridsize1);
    % hosvd
    % [S U sv tol] = hosvd_lpv(data, dep, gridsize, svtol, keep)

    [S1 U1 sv tol] = hosvd_lpv(lpvdata1, dep1, gridsize1, svtol, keep);

    % generating tight polytopic representation
    U1 = genhull(U1, hull);
    S1 = coretensor(U1, lpvdata1, dep1);

    % check model approximation error
    [maxerr meanerr] = tperror(lpvux, S1, U1, domain1, 1000);
    % disp('max and mean error:'); 
    % disp(maxerr);
    % disp(meanerr);
end