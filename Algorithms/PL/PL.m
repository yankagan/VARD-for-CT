% The penalized likelihood algorithm implemented here has not been presented in any
% published paper but the possiblity of using it was mentioned in the paper
%``Alternating Minimization Algorithm for Transmission Tomography" by O'Sullivan and Benac,
% published in IEEE Trans. on Medical Imaging, 2007.
% It is discussed in more detail in ``Alternating Minimization Algorithm
% with Automatic Relevance Determination for Transmission Tomography under Poisson Noise"
% where it is called a Maximum a Posteriori (MAP) estimator
% Code developed by Y. Kaganovsky and S. Degirmenci,  Dec. 2014
% 1D trust region method was suggested by I. Odinaka based on
% the paper ``Newton's Method with a Model Trust Region Modification'' by
% D. C.  Sorensen, published in SIAM Journal of Numerical Analysis, 1982


%%%%%%%%%%%%  parameters for the Newton trust region method used when minimizing the 1D Surrogates  %%%%%%%%%%%%%%%%%%
max_iter = 15;                  % maximal number of Newton steps
step_Tol = 1e-10;               % Threshold for step size
eta1 = 0.1;                     % define ratio = (actual reduction)/(predicted reduction based on quadratic apporximation); when ration<eta1, this is considered to be a bad prediction
eta2 = 0.9;                     % when ratio>eta2 this is considered to be a very good prediction
gam1 = 0.5;                     % the factor by which to shrink the trust region to the left when a bad prediction is identified
gam2 = 2;                       % the factor by which to extend the trust region from the right when a good prediction is identified

x = zeros(N,1);                 % initial guess for posterior mean
obj_fun = zeros(1,N_iter);      % vector to save the history of the objective

del = 1./del;                   % in this file,  delta is defined as 1/delta in the manuscript
reverseStr  =  '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

for iter = 1:N_iter
    msg  =  sprintf('Now running iteration No. %d/%d', iter,N_iter);
    fprintf([reverseStr, msg]);
    reverseStr  =  repmat(sprintf('\b'), 1, length(msg));
    Proj = H_t.'*x;                                                                     % (line 4 of Algorithm 3 in the manuscript), this is faster than H*x due to the way Matlab stores the matrix in memory
    E = eta.*exp(-Proj);                                                                % (line 5 of Algorithm 3 in the manuscript)
    BP_E = H.'*E;                                                                       % (line 6 of Algorithm 3 in the manuscript)
    sv = Psi_v*x;                                                                       % vertical differences
    sh = Psi_h*x;                                                                       % horizontal differences
    obj_fun(iter) = y.'*Proj+sum(E)+(beta./(del.^2)).*sum(abs(sh.*del)-log(1+abs(sh.*del)))+sum(abs(sv.*del)-log(1+abs(sv.*del)));   % evaluate objective at current iteration
    % update estimate using Newton-trust region method (decrease the surrogate function on line 9 of Algorithm 3)%
    % Note that Matlab automatically exploits multiple CPU cores and parallel computing when performing the vectorized operations below
    newton_iter = 0;                                                                    % initialize iteration number
    rad = ones(numel(x),1);                                                             % initialize trust region radius
    fun_hand2  =  @(arg)surrogate_PL(arg,x,y,BP_y,BP_E,eta,Z_ml,Psi_h,Psi_v,beta,del);  % This function computes the 1D surrogate and its derivatives
    [fun der1 der2] =  fun_hand2(x);                                                    % compute the surrogate function and its derivatives at the initial guess
    step = -der1./der2;                                                                 % Newton Step
    ind_keep = find((abs(step)>step_Tol)&(rad>step_Tol)&(newton_iter<max_iter));        % indexes of pixels to be updated
    while ~isempty(ind_keep)
        newton_iter = newton_iter+1;                                                    % increase iteration number
        step(find(abs(step)>rad)) = rad(find(abs(step)>rad));                           % find steps beyond trust region and reduce them to be inside trust region
        new_fun = fun_hand2(x+step);                                                    % evaluate the surrogate at the new point
        act_red = fun-new_fun;                                                          % actual reduction in the surrogate function
        pred_red = -der1.*step-0.5.*der2.*step.^2;                                      % predicted reduction in the surrogate function
        ratio = act_red./pred_red;                                                      % this ratio will determine whether the quadratic approximation is good or not
        % cases when the predicted decrease does'nt agree with the actual decrease
        ind_bad = find((ratio<eta1)|(act_red<0));                                       % find pixels with steps that result in bad prediction
        ind_small = find(abs(step(ind_bad))<rad(ind_bad));                              % within the above find pixels within the trust region boundary
        rad(ind_bad(ind_small)) = abs(step(ind_bad(ind_small)));                        % shrink trust region to the current step
        rad(ind_bad) = rad(ind_bad).*gam1;                                              % shrink the trust region some more
        % cases when the predicted decrease agrees well with the actual decrease
        ind_good = find((ratio >= eta1)&(act_red > 0));                                 % find pixels with good prediction
        ind_update = intersect(ind_keep,ind_good);                                      % which pixels to update
        x(ind_update) = x(ind_update)+step(ind_update);                                 % update the solution
        fun_old = fun;                                                                  % save old value of surrogate
        [fun der1 der2] =  fun_hand2(x);                                                % compute the surrogate function and its derivatives at the new value
        if ~isempty(find(fun>fun_old))                                                  % Inform the user If there are any increases in the 1D surrogates
            fprintf(['One of the 1D surrogates increased at iteration ' num2str(iter) '\nPlease try to change the parametes of the Newton trust region method. \n']);
        end
        step = -der1./der2;                                                             % compute next Newton step
        % cases when the prediction is very good
        ind_very_good = find((ratio>eta2)&(act_red>0));                                 % indexes with very good prediction
        rad(ind_very_good) = rad(ind_very_good).*gam2;                                  % extend the trust region
        % indexes of the pixels that need to be updated
        ind_keep = find((abs(step)>step_Tol)&(rad>step_Tol)&(newton_iter<max_iter));    % update indexes of pixels to be updated in the next iteration
    end
    x(find(x<0)) = 0;                                                                   % Threshold values to be positive
end
time = toc;
fprintf(['\nFinished! total run time was ' num2str(round(time./60)) ' mins \n']);       % display total run time

%%%  plot figures     %%%
% plot reconstructed image %
x_plot = zeros(img_size);
x_plot(ind_circ) = x;
fig = figure;
imagesc(x_plot);
colormap gray
caxis([0.15 0.35]);
colorbar;
title(['Reconstructed Image | PL after ' num2str(N_iter) ' Iterations with \beta = ' num2str(beta) ', \delta = ' num2str(del)],'fontsize', 14);
set(gca, 'fontsize', 14)
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

% plot objective history %
fig = figure;
semilogx(1:iter,obj_fun,'linewidth',3);
title(['Objective vs Iteration Number | PL \beta = ' num2str(beta) ', \delta = ' num2str(del)],'fontsize', 14);
set(gca, 'fontsize', 14)
xlabel('Iteration Number (log-scale)');
ylabel('Objective function value');
xlim([0 N_iter]);
set(gca, 'fontsize', 14)
grid on
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

RMSE_PL = norm(x-object(ind_circ))./norm(object).*100;
fprintf(['The RMSE for PL with beta = ' num2str(beta) ' and delta = ' num2str(del) ' is ' num2str(RMSE_PL) '%% \n']);
