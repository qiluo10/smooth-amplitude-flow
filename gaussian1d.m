clear
close all

Params.n2            = 1;
Params.n1            = 1000;  % signal dimension
Params.no_ini        = 0;     % no initialization?

Params.cplx_flag     = 0;     % set 0 for real signal


if Params.cplx_flag  == 0;
    Params.m             = floor(1.6* Params.n1) ;
else
    Params.m             = floor(2.8* Params.n1) ;
end

Params.T             = 5000;  	% number of gradient iterations
Params.npower_iter   = 250; 	% number of power iterations
Params.alpha         = 0.5;     % weighting parameter in the initialization


Params.mu         = 4* (1 - Params.cplx_flag) + 7* Params.cplx_flag;
Params.order         = 4;
Params.epsilon       = 10e-1;
Params.trun_err      = 1e-5;
Params.SNR = inf;

fprintf('saf: ')
    
    


%% Make signal and data (noiseless)

Amatrix = (1 * randn(Params.m, Params.n1) + Params.cplx_flag * 1i * randn(Params.m, Params.n1)) / (sqrt(2)^Params.cplx_flag);
x       = 1 * randn(Params.n1, 1) + Params.cplx_flag * 1i * randn(Params.n1, 1);

A       = @(I) Amatrix  * I;
At      = @(Y) Amatrix' * Y;

% var     = 0e-1;
noise   = randn(Params.m, 1);
y       = (abs(A(x)));

noipow  = norm(y).^2/norm(noise).^2 / 10^(Params.SNR/10);
y       = (y+ sqrt(noipow)*noise).^2;



%% run SAF algorithm


start_time=tic;

[outs, z] = saf1d(y, x, Params, Amatrix); 

dt =toc(start_time);

Relerrs_RAF = outs.Relerrs;
taus        = outs.taus;
iter_num   = length(Relerrs_RAF);






    
    
    relerrs = outs.Relerrs;
    
    fprintf('ini_err: %.2e, final_err: %.2e\n',...
        relerrs(1),relerrs(end))

% close all

%%

%%
