% The VARD Algorithm implemented here has been presented in the paper
%``Alternating Minimization Algorithm with Automatic Relevance Determination for Transmission Tomography under Poisson Noise
% by Y. Kaganovsky, S. Han, S. Degirmenci,  D. Politte, D. Brady, J. O'Sullivan and L. Carin.

% Version 1.0, Dec. 2014
% Code developed by Y. Kaganovsky and S. Degirmenci
% 1D trust region method was suggested by I. Odinaka based on the paper ``Newton's Method with a Model Trust Region Modification by D. C. Sorensen
% The modification that separates the trust region into left and right regions was suggested by Y. Kaganovsky

%%%%%%%%%%%%%  Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma0 = 1e2.*ones(N,1);        % initialize prior variances
gamma = [ gamma0; gamma0];      % using the neighborhood matrix in Eq. (6.1) of the manuscript
m = zeros(N,1);                 % initial guess for posterior mean
v   = ones(N,1);                % initial guess for posterior variance

%%%%%%%%%%%%  Parameters for the Newton trust region method used when minimizing the 1D Surrogates for the Variances S_{v_j}  %%%%%%%%%%%%%%%%%%
max_iter = 15;                  % maximal number of Newton steps
step_Tol = 1e-10;               % Threshold for step size
eta1 = 0.1;                     % define ratio = (actual reduction)/(predicted reduction based on quadratic apporximation); when ration<eta1, this is considered to be a bad prediction
eta2 = 0.9;                     % when ratio>eta2 this is considered to be a very good prediction
gam1 = 0.5;                     % the factor by which to shrink the trust region to the left when a bad prediction is identified
gam2 = 2;                       % the factor by which to extend the trust region from the right when a good prediction is identified
gam3 = 0.9;                     % a factor used to set the left boundary close to 0 when a negative solution is obtained

obj_fun = zeros(1,N_iter);      % vector to save the history of the objective

reverseStr = '';
tic

for iter = 1:N_iter
    msg = sprintf('Now running iteration No. %d/%d', iter,N_iter);         % message for user regarding the current iteration number
    fprintf([reverseStr, msg]);                                            % see next line
    reverseStr = repmat(sprintf('\b'), 1, length(msg));                    % reverseStr is used to delete the previous message
    g  = (1./gamma).'*abs(Psi);                                            % g^{(t)} in Eq. (3.14) of the manuscript (line 7 of Algorithm 1)
    xi = (Psi.^2).'*(1./gamma);                                            % \xi^{(t)} in Eq. (3.15) of the manuscript (line 7of Algorithm 1)
    proj_m = H_t.'*m;                                                      % This is faster than H*m due to the way Matlab stores the matrix in memory  (line 4 of Algorithm 1 in the manuscript)
    proj_v = H_sq_t.'*v;                                                   % This is faster than H.^2*v due to the same reason (line 5 of Algorithm 1 in the manuscript)
    E = eta.*exp(-proj_m+0.5.*proj_v);                                     % Predicted Poisson rate based on its expectation with respect to the approximate posterior (line 6 of Algorithm 1 in the manuscript)
    BP_E = H.'*E;                                                          % Backprojection of the predicted Poisson rate for the measurements (line 8 of Algorithm 1 in the manuscript)
    SBP_E = H_sq.'*E;                                                      % Squared Backprojection of the predicted Poisson rate (line 9 of Algorithm 1 in the manuscript)
    Psi_m = Psi*m;                                                         % apply difference operator on mean
    f = (Psi_m./gamma).'*Psi;                                              % f^{(t)} in Eq. (3.13) of the manuscript (line 7 of Algorithm 1)
    % Compute the objective
    fun1 = y.'*proj_m+sum(E)+(1./gamma).'*((Psi_m).^2)./2;
    fun2 = -sum(log(v))./2+(1./gamma).'*(Psi.^2*v)./2+sum(log(gamma))./2;
    obj_fun(iter) = fun1+fun2;
    % update mean with one Newton step (decrease the surrogate function S_{m_j} on line 11 of Algorithm 1)%
    %Note that Matlab automatically exploits multiple CPU cores when performing the vectorized operations below
    der = BP_y-BP_E+f.';                                                   % first derivative of the 1D surrogate for the mean S_{m_j}
    der2 = Z_vard.*BP_E+2.*g.';                                            % second derivative of the above
    m = m-der./der2;                                                       % Take one Newton step
    m(find(m<0)) = 0;                                                      % Threshold values to be positive
    % update variances using trust region Newton method (decrease the surrogate function S_{v_j} on line 12 of Algorithm 1) %
    fun_hand = @(arg)surrogate_VARD(arg,v,SBP_E,xi,Z_vard);                % This function computes the 1D surrogate for the variance S_{v_j} as well as its derivatives
    newton_iter = 0;                                                       % initialize iteration number
    rad_r = ones(numel(v),1);                                              % set right boundary of trust region
    rad_l = -0.9.*v;                                                       % set left boundary of trust region to exclude negatives
    [fun der1 der2] = fun_hand(v);                                         % evaluate the surrogate and its derivatives at the current guess
    step = -der1./der2;                                                    % Newton step
    ind_keep = find((abs(step)>step_Tol)&((rad_r>step_Tol)|(abs(rad_l)>step_Tol))&(newton_iter<max_iter));            % indexes of pixels to be updated
    % Newton-trust region %
    %Note that Matlab automatically exploits multiple CPU cores and parallel computing when performing the vectorized operations below
    while ~isempty(ind_keep)
        newton_iter = newton_iter+1;                                                                                  % increase iteration number
        step(find(step>rad_r)) = rad_r(find(step>rad_r));                                                             % find steps beyond right boundary and reduce them to be inside the trust region
        step(find(step<rad_l)) = rad_l(find(step<rad_l));                                                             % find steps beyond left boundary and reduce them to be inside the trust region
        new_fun = fun_hand(v+step);                                                                                   % evaluate the surrogate at the new point
        act_red = fun-new_fun;                                                                                        % actual reduction in the surrogate function
        pred_red = -der1.*step-0.5.*der2.*step.^2;                                                                    % predicted reduction in the surrogate function
        ratio = act_red./pred_red;                                                                                    % this ratio will determine whether the quadratic approximation is good or not
        % cases when the predicted decrease  does'nt agree with the actual decrease
        ind_bad_right = find(((ratio<eta1)|(act_red<0))&(step>0));                                                    % find pixels with steps to the right that result in bad prediction
        ind_small = find(step(ind_bad_right)<rad_r(ind_bad_right));                                                   % within the above, find pixels not going beyond the right boundary
        rad_r(ind_bad_right(ind_small)) = step(ind_bad_right(ind_small));                                             % shrink trust region from the right to the current step
        rad_r(ind_bad_right) = rad_r(ind_bad_right).*gam1;                                                            % shrink the trust region from the right some more
        ind_bad_left = find(((ratio<eta1)|(act_red<0))&(step<= 0));                                                  % find pixels with steps to the left that result in bad prediction
        ind_out = find((step(ind_bad_left)<rad_l(ind_bad_left))|((v(ind_bad_left)+step(ind_bad_left))<0));            % within the above, find pixels with negative values or going beyond the left boundary
        rad_l(ind_bad_left(ind_out)) = max([step(ind_bad_left(ind_out)) -v(ind_bad_left(ind_out)).*gam3],[],2);       % shrink the trust region from the  left to the current step or close to zero if the step results in a negative solution
        rad_l(ind_bad_left) = rad_l(ind_bad_left).*gam1;                                                              % shrink the trust region from the left
        % cases when the predicted decrease agrees well with the actual decrease
        ind_good = find((ratio>= eta1)&(act_red>= 0));                                                              % find pixels with good prediction
        ind_neg = find((v(ind_good)+step(ind_good))<0);                                                               % within the above, find the ones with negative values (should not happen, but just to be safe)
        rad_l(ind_good(ind_neg)) = -v(ind_good(ind_neg)).*gam3;                                                       % set the left boundary close to 0
        ind_pos = find((v(ind_good)+step(ind_good))>= 0);                                                            % find non-negative values within pixels with good prediction
        ind_update = intersect(ind_keep,ind_good(ind_pos));                                                           % indexes of pixels to be update
        v(ind_update) = v(ind_update)+step(ind_update);                                                               % update variances
        fun_prev = fun;                                                                                               % save last solution
        [fun der1 der2] = fun_hand(v);                                                                                % evaluate surrogate at the current solution
        if ~isempty(find(fun>fun_prev))                                                                               % Inform the user if there are any increases in the 1D surrogates
            fprintf(['One of the 1D surrogates for the variances increased at iteration ' num2str(iter) '\nPlease try to change the parametes of the Newton trust region method. \n']);
        end
        step = -der1./der2;                                                                                           % Newton step
        % cases when the prediction is very good
        ind_very_good = find((ratio>eta2)&(step>0)&(act_red>0));                                                      % indexes with very good prediction
        rad_r(ind_very_good) = rad_r(ind_very_good).*gam2;                                                            % extend the trust region from the right
        % indexes of the pixels that need to be updated
        ind_keep = find((abs(step)>step_Tol)&((rad_r>step_Tol)|(abs(rad_l)>step_Tol))&(newton_iter<max_iter));        % update indexes of pixels to be updated in the next iteration
    end
    %%%%%%%%   update gamma   %%%%%%%%%%%%
    gamma0 = (Psi_h*m).^2+Psi_h.^2*v+(Psi_v*m).^2+Psi_v.^2*v ;                                                        % update prior hyperparameters (line 13 of Algorithm 1 in the manuscript)
    gamma = [gamma0; gamma0];
end
time = toc;
fprintf(['\nFinished! total run time was ' num2str(round(time./60)) ' mins \n']);

%%%  plot figures     %%%
m_plot = zeros(img_size);
m_plot(ind_circ) = m;
v_plot = zeros(img_size);
v_plot(ind_circ) = v;
gamma_plot = zeros(img_size);
gamma_plot(ind_circ) = gamma0;

% plot posterior mean %
fig = figure;
imagesc(m_plot);
colormap gray
title(['Reconstructed Image | VARD after ' num2str(N_iter) ' Iterations'],'fontsize', 14);
caxis([0.15 0.35]);
colorbar
set(gca, 'fontsize', 14)
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

% plot posterior std %
fig = figure;
imagesc(sqrt(v_plot));
colormap gray
title(['Posterior Standard Deviation | VARD after ' num2str(N_iter) ' Iterations'],'fontsize', 14);
caxis([0 0.01]);
colorbar
set(gca, 'fontsize', 14)
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

% plot learned hyperparameters %
fig = figure;
imagesc(sqrt(gamma_plot));
colormap gray
title(['Learned Hyperparameters | VARD after ' num2str(N_iter) ' Iterations'],'fontsize', 14);
caxis([0 0.05]);
colorbar
set(gca, 'fontsize', 14)
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

% plot objective history %
fig = figure;
semilogx(1:iter,obj_fun,'linewidth',3);
title('Objective vs Iteration Number | VARD','fontsize', 14);
xlabel('Iteration Number (log-scale)');
ylabel('Objective function value');
set(gca, 'fontsize', 14)
xlim([0 N_iter]);
grid on
set(fig, 'units', 'inches', 'position', [4 1 8 7]);

RMSE_vard = norm(m-object(ind_circ))./norm(object).*100;
fprintf(['The RMSE for VARD is ' num2str(RMSE_vard) '%% \n']);

