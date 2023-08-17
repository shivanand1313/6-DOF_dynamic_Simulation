function dstates = six_dof_model(t, states)

%% params
m = 750;g = 9.81; Ixx = 873; Iyy = 907; Izz = 1680; Ixz = 1144;S = 12.47; b = 10.47; cb = 1.211; 
c_d0 = 0.035; k = 0.045; cL_0 = 0.370; cL_a = 5.0; cL_q = 37.211; cL_de = 0.374;
cm_0 = 0.091; cm_a = -2.937; cm_q = -8.719; cm_de = -0.735;
cy_0 = 0; cy_b = -0.531; cy_p = -0.0571; cy_r = 0.4657; cy_dr = 0.1502;
cl_0 = 0; cl_b = -0.031; cl_p = -0.262; cl_r = -0.0541; cl_dr = 0.005; cl_da = -0.153;
cn_0 = 0; cn_b = 0.01; cn_p = -0.007; cn_da = 0; cn_r = -0.067; cn_dr = -0.047;
%%

statesCell = num2cell(states);
[u, v, w, x, y, z, p, q, r, phi, theta, psi] = statesCell{:};

%% Controlled Input function
ctrl = cont(t, states);

%% AeroF&M__Kinematics 
dele = ctrl(1); dela = ctrl(2); delr = ctrl(3);

rho0 = 1.225; % kg/m^3
H0 = 10400; %m (scale height)

rho = rho0*exp(-z/H0);

V_inf = sqrt(u^2 + v^2 + w^2);
alpha = atan2(w,u);
beta = asin(v/V_inf);

ph = p*b/2/V_inf;    qh = q*cb/2/V_inf;   rh = r*b/2/V_inf;

cL = cL_0 + cL_a*alpha + cL_q*qh + cL_de*dele;
cD = c_d0 + k*cL^2;
cY = cy_0 + cy_b*beta + cy_dr*delr;
cM = cm_0 + cm_a*alpha + cm_q*qh + cm_de*dele;
cl = cl_0 + cl_b*beta + cl_p*ph + cl_r*rh + cl_da*dela + cl_dr*delr;
cN = cn_0 + cn_b*beta + cn_p*ph + cn_r*rh + cn_da*dela + cn_dr*delr;

cX = cL*sin(alpha) - cD*cos(alpha);
cZ = -(cL*cos(alpha) + cD*sin(alpha));

Q = 1/2*rho*V_inf^2;

Fx= Q*S*cX;     Fy= Q*S*cY;     Fz= Q*S*cZ;
l = Q*S*b*cl;   M = Q*S*cb*cM;  N = Q*S*b*cN;

T = ctrl(4);

du = Fx/m + T/m - q*w + r*v - g*sin(theta);
dv = Fy/m - r*u + p*w + g*cos(theta)*sin(phi);
dw = Fz/m - p*v + q*u + g*cos(theta)*cos(phi);

%% Rotation Matrix
Rpsi = [cos(psi) sin(psi) 0;
       -sin(psi) cos(psi) 0;
            0       0     1];

Rtheta = [cos(theta)  0  -sin(theta);
             0        1       0; 
          sin(theta)  0  cos(theta)];

Rphi = [1      0         0;
        0   cos(phi)  sin(phi);
        0  -sin(phi)  cos(phi)];

b2i_mat = Rpsi'*Rtheta'*Rphi';

pqr2eul = [1  sin(phi)*tan(theta)  cos(phi)*tan(theta);
           0      cos(phi)             -sin(phi);
           0  sin(phi)/cos(theta)  cos(phi)/cos(theta)];
%% 
dxx = b2i_mat*[u; v; w];

dx = dxx(1);    dy = dxx(2);    dz = dxx(3);

dp = (l*Izz + N*Ixz - p*q*(Ixz*(Iyy-Ixx-Izz)) - q*r*(Izz*(Izz-Iyy)+Ixz^2))/(Ixx*Izz - Ixz^2);
dq = (M + p*r*(Izz-Ixx) - Ixz*(p^2-r^2))/Iyy;
dr = (l*Ixz + N*Ixx + p*q*(Ixz^2+Ixx*(Ixx-Iyy)) + q*r*(Ixz*(Iyy-Ixx-Izz)))/(Ixx*Izz - Ixz^2);

deul = pqr2eul*[p;q;r];
dphi = deul(1); dtheta = deul(2);   dpsi = deul(3);

dstates = [du, dv, dw, dx, dy, dz, dp, dq, dr, dphi, dtheta, dpsi]';
end