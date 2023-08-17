% prints the trim conditions

%% params
m = 750;g = 9.81; Ixx = 873; Iyy = 907; Izz = 1680; Ixz = 1144;S = 12.47; b = 10.47; cb = 1.211; 
c_d0 = 0.035; k = 0.045; cL_0 = 0.370; cL_a = 5.0; cL_q = 37.211; cL_de = 0.374;
cm_0 = 0.091; cm_a = -2.937; cm_q = -8.719; cm_de = -0.735;
cy_0 = 0; cy_b = -0.531; cy_p = -0.0571; cy_r = 0.4657; cy_dr = 0.1502;
cl_0 = 0; cl_b = -0.031; cl_p = -0.262; cl_r = -0.0541; cl_dr = 0.005; cl_da = -0.153;
cn_0 = 0; cn_b = 0.01; cn_p = -0.007; cn_da = 0; cn_r = -0.067; cn_dr = -0.047;
%%
z = 3000;% m
rho0 = 1.225; % kg/m^3
H0 = 10400; % m (scale height)

rho = rho0*exp(-z/H0);

Vinf = 55.57; % m/s---------------------

W = m*g;
Q = 1/2*rho*Vinf^2;

cLtrim = W / (Q*S)
cDtrim = c_d0 + k*cLtrim^2;
thrust = Q*S*cDtrim

res = ([cL_a cL_de; cm_a cm_de])\[cLtrim-cL_0; -cm_0];

alpha_trim = res(1)
dele_trim = res(2)

u_trim = Vinf*cos(alpha_trim)
w_trim = Vinf*sin(alpha_trim)