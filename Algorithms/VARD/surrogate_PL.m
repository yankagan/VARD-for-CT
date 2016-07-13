function [fun der1 der2] =  surrogate_PL(x,x_prev,y,BP_y,BP_E,y_r,Z,Dh,Dv,beta,del)

% Compute the surrogate for the Huber penalty with its 1st and 2nd derivates %
% The Huber penalty is of the form: beta*sum_{j,k, in neighborhood} P(x_j-x_k)
% where P(x) = (|del*x|-log(1+|del*x|))/(del^2)

%% Inputs %%
% x      -   current estimate of the image
% x_prev -   previous estimate of the image
% y      -   data with object
% BP_y   -   backprojected data
% BP_E   -   backprojection of predicted measurements mean (Poisson rate) 
% y_r    -   data for air scan
% Dh     -   horizontal difference operator (1's on main diagonal and -1 on the diagonal above the main diagonal )
% Dv     -   vertical difference operator (1's on main diagonal and -1 on the Nth diagonal below ; N  = num. of pixels along each image dimension). 
% beta   -   weight of the Huber penalty 
% del    -   scaling constant in the Huber penalty  

%compute horizontal/vertical differences 
Dhm = Dh*x_prev;
Dvm = Dv*x_prev;
Dhtm = Dh.'*x_prev;
Dvtm = Dv.'*x_prev;
arg_h1 = del.*(Dhm+2.*(x-x_prev));
arg_v1 = del.*(Dvm+2.*(x-x_prev));
arg_h2 = del.*(Dhtm+2.*(x-x_prev));
arg_v2 = del.*(Dvtm+2.*(x-x_prev));

% compute surrogate for each pixel
fun_lik = BP_y.*x+BP_E.*exp(-Z.*(x-x_prev))./Z;                                                            % data-fit
fun_pen1 = beta./(2.*del.^2).*(abs(arg_h1)-log(1+abs(arg_h1))+abs(arg_v1)-log(1+abs(arg_v1))); 
fun_pen2 = beta./(2.*del.^2).*(abs(arg_h2)-log(1+abs(arg_h2))+abs(arg_v2)-log(1+abs(arg_v2))); 
fun_pen = fun_pen1+fun_pen2;                                                                               % Huber penalty 
fun = fun_lik+fun_pen;


% compute 1st derivative of surrogate for each pixel
der1_lik = BP_y-BP_E.*exp(-Z.*(x-x_prev));                                                                 % data-fit 
der1_pen1 = (beta./del).*( sign(arg_h1).*(1-1./(1+abs(arg_h1))) + sign(arg_v1).*(1-1./(1+abs(arg_v1))) ); 
der1_pen2 = (beta./del).*( sign(arg_h2).*(1-1./(1+abs(arg_h2))) + sign(arg_v2).*(1-1./(1+abs(arg_v2))) ); 
der1_pen = der1_pen1+der1_pen2;
der1 = der1_lik+der1_pen;


% compute 2nd derivative of surrogate for each pixel
der2_lik = Z.*BP_E.*exp(-Z.*(x-x_prev));                                                                   % data-fit 
der2_pen1 = beta.*( 2./((1+abs(arg_h1)).^2)+ 2./((1+abs(arg_v1)).^2)  );
der2_pen2 = beta.*( 2./((1+abs(arg_h2)).^2)+ 2./((1+abs(arg_v2)).^2)  );
der2_pen = der2_pen1+der2_pen2;                                                                            % Huber penalty 
der2 = der2_lik+der2_pen;


         
