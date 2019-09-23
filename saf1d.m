%% Implementation of the Smooth Amplitude Flow algorithm proposed in the paper
%  Phase Retrieval via Smooth Amplitude Flow
%  by  Q. Luo, H. Wang, and S. Lin from NUDT, China
%  The code below is adapted from implementation of the Wirtinger Flow
%  algorithm implemented by E. Candes, X. Li, and M. Soltanolkotabi.

function [outs, z] = saf1d(y, x, Params, Amatrix)



%% Initialization
if Params.cplx_flag
    z0 = randn(Params.n1, Params.n2) + 1i*randn(Params.n1, Params.n2);
else
    z0 = randn(Params.n1, Params.n2);
end

z      = z0 / norm(z0, 'fro');    % Initial guess
AmatT  = Amatrix';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normest = sqrt(sum(y(:)) / numel(y(:)));    % Estimate norm to scale eigenvector
m       = Params.m;
ymag    = sqrt(y);
order   = Params.order;


epsilon = Params.epsilon;
ymag_order = epsilon^order * ymag.^order;
trun_err= Params.trun_err;


%% The weighted maximal correlation initialization can be computed using power iterations
%% or the Lanczos algorithm, and the latter performs well when m/n is small

if ~Params.no_ini
    [ysort, ~] = sort(y, 'ascend');
    ythresh    = ysort(round(m / 1.3));
    ind        = y >= ythresh;
    Aselect    = Amatrix(ind, :);
    weights    = (ymag(ind)).^(Params.alpha); % weights w_i
    
    for i = 1:Params.npower_iter                   % Power iterations
        z  = Aselect' * (weights .* (Aselect * z));
        z  = z / norm(z, 'fro');
    end
end

z = normest * z;


Relerrs = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro'); % Initial rel. error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tt = (1+epsilon^order)^(1/order)*ymag;


[f,gradf] = SAF_Objective; 


beta = 0.2;


for t = 1: Params.T
    tau = Params.mu;
    
    
    Az    = Amatrix * z;
    f_now = f(Az);
    grad_now = AmatT*gradf(Az);
    
    count = 0;
    ngrad2 = norm(grad_now)^2;
    
    while count<=1
        znew  = z - tau  * grad_now;
        Aznew = Amatrix * znew;
        f_znew = f(Aznew);
        
        
        if f_znew-f_now <= - tau*beta*ngrad2
            
            %             disp(count)
            break;
        end
        
        
        
        tau = tau*0.2;
        count = count+1;
    end
    z = znew;
    %     disp(tau)
    %     if count==0
    %         tau = tau*1.4;
    %     end
    
    %     if tau<2e-15
    %         break
    %     end
    
    errcnt = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro');
    
    Relerrs = [Relerrs; errcnt]; %#ok<AGROW>
    taus(t)    = tau;
    
    
    if errcnt<= trun_err
        break
    end
end


outs.Relerrs = Relerrs;
outs.taus    = taus; % output the stepth at each iteration



    function [f, gradf] = SAF_Objective(~, ~)
        f = @(z) 0.5 * norm( (abs(z).^order+ymag_order).^(1/order)-tt   ).^2/m;
        gradf = @(z) ((abs(z).^order+ymag_order).^(1/order)-tt).*(abs(z).^order+...
            ymag_order).^(1/order-1).* abs(z).^(order-1).*sign(z)/m;
    end






end




