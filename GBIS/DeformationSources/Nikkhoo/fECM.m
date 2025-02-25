function [ue,un,uv,dV,DV,Ns]=fECM(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,...
    ax,ay,az,p,mu,lambda,DepthRef,Nmax,Cr)
% fECM
% calculates the surface displacements caused by a uniformly-pressurized 
% finite ellipsoidal cavity in a uniform elastic half-space.
%
% fECM: finite Ellipsoidal Cavity Model
% pCDM: point Compound Dislocation Model
% PTD: Point Tensile Dislocation
% EFCS: Earth-Fixed Coordinate System
%
% INPUTS
% X and Y:
% Horizontal coordinates of calculation points in EFCS (East, North, Up).
% The matrices X and Y must have the same size.
%
% X0 and Y0:
% Horizontal coordinates (in EFCS) of the center of the cavity. X0 and Y0 
% have the same unit as X and Y.
%
% depth:
% The depth to the "reference point", that is, the center or top point of
% the cavity, specified by the "DepthRef" input. The depth is a positive
% number with the same unit as X0 and Y0.
% 
% omegaX, omegaY and omegaZ:
% Clockwise rotation angles about the X, Y and Z axes, respectively, which
% specify the orientation of the finite ECM in space. These angles must be
% given in degrees.
%
% ax, ay and az:
% Semi-axes of the finite ECM along the X, Y and Z axes, respectively,
% before applying the rotations. ax, ay and az have the same unit as X and
% Y.
%
% p:
% Pressure on the cavity walls. p has the same unit as the Lamé constants.
%
% mu and lambda:
% Lamé constants.
%
% DepthRef:
% is either 'C' or 'T' and specifies the "reference point" for the "depth".
% If 'C' is chosen, depth = dC (depth to the center), and if 'T' is chosen,
% depth = dT (depth to the cavity top point).
%
% OPTIONAL INPUTS
% Nmax:
% Maximum total number of allowed point CDMs. If not given as an input the
% default value Nmax = 5e3 is used.
%
% Cr:
% The "grid spacing parameter" (see reference paper 1). If not given as an
% input the default value Cr = 14 is used.
%
%
% OUTPUTS
% ue, un and uv:
% Calculated displacement vector components in EFCS. ue, un and uv have the
% same unit as X, Y and Z.
%
% dV and DV:
% Volume change and potency of the finite ECM. dV and DV have the same
% unit, that is, the unit of all the other quantities with dimension
% distance (i.e., X,Y,depth,ax,ay,az,ue,un,uv) to the power of 3.
%
% Ns:
% The total number of point CDMs within the cavity.
%
%
% Example: Calculate and plot surface displacements along the X axis
%
% [X,Y] = meshgrid(0:50:7000,0);
% X0 = 0; Y0 = 0; dT = 550; omegaX = 11; omegaY = -7; omegaZ = 120;
% ax = 450; ay = 600; az = 1225; p = 1e6; mu = 1e9; lambda = 1e9;
% [ue,~,uv,dV,DV,Ns] = fECM(X,Y,X0,Y0,dT,omegaX,omegaY,omegaZ,ax,ay,az,...
%     p,mu,lambda,'T',4000,12);
% figure
% plot(X,uv*1e2,'k')
% hold on
% plot(X,ue*1e2,'b')
% title('Displacements')
% xlabel('X (meter)')
% ylabel('u_{i} (cm)')
% legend('u_v','u_e')

% Reference journal articles:
% 1)
% Nikkhoo, M., Rivalta, E. (2023):
% Surface deformations and gravity changes caused by pressurized finite 
% ellipsoidal cavities. Geophysical Journal International, doi: 
% 10.1093/gji/ggac351
%
% 2)
% Nikkhoo, M., Rivalta, E. (2022):
% Analytical solutions for gravity changes caused by triaxial volumetric
% sources. Geophysical Research Letters, doi: 10.1029/2021GL095442
%
% 3)
% Nikkhoo, M., Walter, T. R., Lundgren, P. R., Prats-Iraola, P. (2017):
% Compound dislocation models (CDMs) for volcano deformation analyses.
% Geophysical Journal International, doi: 10.1093/gji/ggw427
%
% 4)
% Eshelby, J. D. (1957):
% The determination of the elastic field of an ellipsoidal inclusion, and
% related problems.
% Proceedings of the royal society of London. Series A. Mathematical and
% physical sciences. 241 (1226), 376-396. doi: 10.1098/rspa.1957.0133
%
% 5)
% Carlson, B. C. (1995):
% Numerical computation of real or complex elliptic integrals.
% Numer. Algor., 10(1), 13–26. doi: 10.1007/BF02198293
%
% 6)
% Klein, P. P. (2012):
% On the ellipsoid and plane intersection equation.
% Applied Mathematics, 3 (11), 1634-1640. doi: 10.4236/am.2012.311226

% Copyright (c) 2022 Mehdi Nikkhoo
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction, including
% without limitation the rights to use, copy, modify, merge, publish,
% distribute, sublicense, and/or sell copies of the Software, and to permit
% persons to whom the Software is furnished to do so, subject to the
% following conditions:
%
% The above copyright notice and this permission notice shall be included
% in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
% NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
% OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
% USE OR OTHER DEALINGS IN THE SOFTWARE.

% I appreciate any comments or bug reports.

% Mehdi Nikkhoo
% Created: 2021.7.24
% Last modified: 2022.10.24
%
% Section 2.1, Physics of Earthquakes and Volcanoes
% Department 2, Geophysics
% Helmholtz Centre Potsdam
% German Research Centre for Geosciences (GFZ)
%
% email:
% mehdi.nikkhoo@gfz-potsdam.de
% mehdi.nikkhoo@gmail.com
%
% website:
% http://www.volcanodeformation.com

if nargin<16
    Nmax = 5e3;
    Cr = 14;
elseif nargin<17
    Cr = 14;
end

[Zt,Zb] = EllTopBot(depth,omegaX,omegaY,omegaZ,ax,ay,az);
if strcmpi(DepthRef,'T')
    aV = (Zt-Zb)/2;
    depth = depth+aV; % Note that dC = dT+aV!
    Ztop = Zt-aV;
elseif strcmpi(DepthRef,'C')
    Ztop = Zt;
else
    error('Undefined input argument for ''RefPoint''!')
end

if Ztop>=0
    error('Input error: the cavity is too shallow!')
end

X = X(:);
Y = Y(:);

nu = lambda/(lambda+mu)/2; % Poisson's ratio
K = lambda+2*mu/3; % Bulk modulus

r0 = 1e-12; % Threshold for stable shape tensor calculation
ax(ax<r0) = r0;
ay(ay<r0) = r0;
az(az<r0) = r0;
[ai,Ind] = sort([ax ay az],'descend');

S = ShapeTensorECM(ai(1),ai(2),ai(3),nu);
Sm = [S(1)-1 S(2) S(3);S(4) S(5)-1 S(6);S(7) S(8) S(9)-1]; % Shape tensor

eT = -inv(Sm)*p*ones(3,1)/3/K; % Transformation strain
% For uniformly-pressurized ellipsoids the eT elements have the same sign!
eT(sign(eT)~=sign(p)) = 0;
V = 4/3*pi*ax*ay*az; % Ellipsoid volume

% calculate actual volume change
dV = (sum(eT)-p/K)*V;
% calculate potency
DV = dV+p*V/K; % Also DV = sum(M)/3/K could be used!

tmp = zeros(3,1);
tmp(Ind) = V*eT;
DVx = tmp(1);
DVy = tmp(2);
DVz = tmp(3);

[Xs,Ys,Zs,Ws] = AdaptiveGrid(X0,Y0,depth,omegaX,omegaY,omegaZ,ax,ay,az,...
    Cr,Nmax);
Ns = numel(Xs);
DVx = DVx*Ws;
DVy = DVy*Ws;
DVz = DVz*Ws;

[ue,un,uv] = pCDM_Vec(X,Y,Xs,Ys,-Zs,omegaX,omegaY,omegaZ,DVx,DVy,DVz,nu);


function [ue,un,uv]=pCDM_Vec(X,Y,X0,Y0,depth,omegaX,omegaY,omegaZ,...
    DVx,DVy,DVz,nu)
% pCDM_Vec
% is the vectorized version of the pCDM function (see the codes published 
% together with reference #3). 

X = X(:);
Y = Y(:);

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

Vstrike1 = [-R(2,1),R(1,1),0];
Vstrike1 = Vstrike1/norm(Vstrike1);
strike1 = atan2(Vstrike1(1),Vstrike1(2))*180/pi;
if isnan(strike1)
    strike1 = 0;
end
dip1 = acosd(R(3,1));

Vstrike2 = [-R(2,2),R(1,2),0];
Vstrike2 = Vstrike2/norm(Vstrike2);
strike2 = atan2(Vstrike2(1),Vstrike2(2))*180/pi;
if isnan(strike2)
    strike2 = 0;
end
dip2 = acosd(R(3,2));

Vstrike3 = [-R(2,3),R(1,3),0];
Vstrike3 = Vstrike3/norm(Vstrike3);
strike3 = atan2(Vstrike3(1),Vstrike3(2))*180/pi;
if isnan(strike3)
    strike3 = 0;
end
dip3 = acosd(R(3,3));

% Calculate contribution of the first PTD
[ue1,un1,uv1] = PTD_disp_Surf_Vec(X,Y,X0,Y0,depth,strike1,dip1,DVx,nu);

% Calculate contribution of the second PTD
[ue2,un2,uv2] = PTD_disp_Surf_Vec(X,Y,X0,Y0,depth,strike2,dip2,DVy,nu);

% Calculate contribution of the third PTD
[ue3,un3,uv3] = PTD_disp_Surf_Vec(X,Y,X0,Y0,depth,strike3,dip3,DVz,nu);

ue = ue1+ue2+ue3;
un = un1+un2+un3;
uv = uv1+uv2+uv3;

function [ue,un,uv]=PTD_disp_Surf_Vec(X,Y,X0,Y0,depth,strike,dip,DV,nu)
% PTD_disp_Surf_Vec
% is the vectorized version of PTD_disp_Surf function in the pCDM function
% (see the codes published together with reference #3).

X = X(:);
Y = Y(:);
% X0, Y0, depth and DV must be matrices of the same size
X0 = X0(:)';
Y0 = Y0(:)';
depth = depth(:)';
DV = DV(:)';

xM = repmat(X,1,numel(X0))-repmat(X0,numel(X),1); % a matrix!
yM = repmat(Y,1,numel(X0))-repmat(Y0,numel(X),1);
d = repmat(depth,numel(X),1);
DV = repmat(DV,numel(X),1);

beta = strike-90;
x = xM*cosd(beta)-yM*sind(beta);
y = xM*sind(beta)+yM*cosd(beta);

r = (x.^2+y.^2+d.^2).^0.5;
q = y*sind(dip)-d*cosd(dip);

I1 = (1-2*nu)*y.*(1./r./(r+d).^2-x.^2.*(3*r+d)./r.^3./(r+d).^3);
I2 = (1-2*nu)*x.*(1./r./(r+d).^2-y.^2.*(3*r+d)./r.^3./(r+d).^3);
I3 = (1-2*nu)*x./r.^3-I2;
I5 = (1-2*nu)*(1./r./(r+d)-x.^2.*(2*r+d)./r.^3./(r+d).^2);

% Note: For a PTD M0 = DV*mu!
ue = DV/2/pi.*(3*x.*q.^2./r.^5-I3*sind(dip)^2);
un = DV/2/pi.*(3*y.*q.^2./r.^5-I1*sind(dip)^2);
uv = DV/2/pi.*(3*d.*q.^2./r.^5-I5*sind(dip)^2);

ue0 = sum(ue,2);
un0 = sum(un,2);
uv = sum(uv,2);

ue = ue0*cosd(beta)+un0*sind(beta);
un = -ue0*sind(beta)+un0*cosd(beta);


function [S]=ShapeTensorECM(a1,a2,a3,nu)
% ShapeTensorECM
% calculates the Eshelby (1957) shape tensor components.

if all([a1 a2 a3]==0)
    S = zeros(9,1);
    return
end

% Calculate Ik and Iij terms for triaxial, oblate and prolate ellipsoids
if a1>a2 && a2>a3 && a3>0
    % General case: triaxial ellipsoid
    sin_theta = sqrt(1-a3^2/a1^2);
    k = sqrt((a1^2-a2^2)/(a1^2-a3^2));
    
    % % Calculate Legendre's incomplete elliptic integrals of the first and
    % % second kind using MATLAB Symbolic Math Toolbox
    % F =  mfun('EllipticF',sin_theta,k);
    % E =  mfun('EllipticE',sin_theta,k);
    
    % Calculate Legendre's incomplete elliptic integrals of the first and
    % second kind using Carlson (1995) method (see Numerical computation of
    % real or complex elliptic integrals. Carlson, B.C. Numerical
    % Algorithms (1995) 10: 13. doi:10.1007/BF02198293)
    tol = 1e-16;
    c = 1/sin_theta^2;
    F = RF(c-1,c-k^2,c,tol);
    E = F-k^2/3*RD(c-1,c-k^2,c,tol);
    
    I1 = 4*pi*a1*a2*a3/(a1^2-a2^2)/sqrt(a1^2-a3^2)*(F-E);
    I3 = 4*pi*a1*a2*a3/(a2^2-a3^2)/sqrt(a1^2-a3^2)*...
        (a2*sqrt(a1^2-a3^2)/a1/a3-E);
    I2 = 4*pi-I1-I3;
    
    I12 = (I2-I1)/(a1^2-a2^2);
    I13 = (I3-I1)/(a1^2-a3^2);
    I11 = (4*pi/a1^2-I12-I13)/3;
    
    I23 = (I3-I2)/(a2^2-a3^2);
    I21 = I12;
    I22 = (4*pi/a2^2-I23-I21)/3;
    
    I31 = I13;
    I32 = I23;
    I33 = (4*pi/a3^2-I31-I32)/3;
    
elseif a1==a2 && a2>a3 && a3>0
    % Special case-1: Oblate ellipsoid
    I1 = 2*pi*a1*a2*a3/(a1^2-a3^2)^1.5*(acos(a3/a1)-a3/a1*...
        sqrt(1-a3^2/a1^2));
    I2 = I1;
    I3 = 4*pi-2*I1;
    
    I13 = (I3-I1)/(a1^2-a3^2);
    I11 = pi/a1^2-I13/4;
    I12 = I11;
    
    I23 = I13;
    I22 = pi/a2^2-I23/4;
    I21 = I12;
    
    I31 = I13;
    I32 = I23;
    I33 = (4*pi/a3^2-2*I31)/3;
    
elseif a1>a2 && a2==a3 && a3>0
    % Special case-2: Prolate ellipsoid
    I2 = 2*pi*a1*a2*a3/(a1^2-a3^2)^1.5*(a1/a3*sqrt(a1^2/a3^2-1)-...
        acosh(a1/a3));
    I3 = I2;
    I1 = 4*pi-2*I2;
    
    I12 = (I2-I1)/(a1^2-a2^2);
    I13 = I12;
    I11 = (4*pi/a1^2-2*I12)/3;
    
    I21 = I12;
    I22 = pi/a2^2-I21/4;
    I23 = I22;
    
    I32 = I23;
    I31 = I13;
    I33 = (4*pi/a3^2-I31-I32)/3;
end

% Calculate the shape-tensor components
if a1==a2 && a2==a3
    % Special case-3: Sphere
    S1111 = (7-5*nu)/15/(1-nu);
    S1122 = (5*nu-1)/15/(1-nu);
    S1133 = (5*nu-1)/15/(1-nu);
    S2211 = (5*nu-1)/15/(1-nu);
    S2222 = (7-5*nu)/15/(1-nu);
    S2233 = (5*nu-1)/15/(1-nu);
    S3311 = (5*nu-1)/15/(1-nu);
    S3322 = (5*nu-1)/15/(1-nu);
    S3333 = (7-5*nu)/15/(1-nu);
else
    % General triaxial, oblate and prolate ellipsoids
    S1111 = 3/8/pi/(1-nu)*a1^2*I11+(1-2*nu)/8/pi/(1-nu)*I1;
    S1122 = 1/8/pi/(1-nu)*a2^2*I12-(1-2*nu)/8/pi/(1-nu)*I1;
    S1133 = 1/8/pi/(1-nu)*a3^2*I13-(1-2*nu)/8/pi/(1-nu)*I1;
    S2211 = 1/8/pi/(1-nu)*a1^2*I21-(1-2*nu)/8/pi/(1-nu)*I2;
    S2222 = 3/8/pi/(1-nu)*a2^2*I22+(1-2*nu)/8/pi/(1-nu)*I2;
    S2233 = 1/8/pi/(1-nu)*a3^2*I23-(1-2*nu)/8/pi/(1-nu)*I2;
    S3311 = 1/8/pi/(1-nu)*a1^2*I31-(1-2*nu)/8/pi/(1-nu)*I3;
    S3322 = 1/8/pi/(1-nu)*a2^2*I32-(1-2*nu)/8/pi/(1-nu)*I3;
    S3333 = 3/8/pi/(1-nu)*a3^2*I33+(1-2*nu)/8/pi/(1-nu)*I3;
end

S = [S1111 S1122 S1133 S2211 S2222 S2233 S3311 S3322 S3333]';

function [rf]=RF(x,y,z,r)
% RF
% calculates the RF term in the Carlson (1995) method for calculating
% elliptic integrals

% r = 1e-16

if any([x,y,z]<0)
    error('x, y and z values must be positive!')
elseif nnz([x,y,z])<2
    error('At most one of the x, y and z values can be zero!')
end

xm = x;
ym = y;
zm = z;
A0 = (x+y+z)/3;
Q = max([abs(A0-x),abs(A0-y),abs(A0-z)])/(3*r)^(1/6);
n = 0;
Am = A0;
while abs(Am)<=Q/(4^n)
    lambdam = sqrt(xm*ym)+sqrt(xm*zm)+sqrt(ym*zm);
    Am = (Am+lambdam)/4;
    xm = (xm+lambdam)/4;
    ym = (ym+lambdam)/4;
    zm = (zm+lambdam)/4;
    n = n+1;
end
X = (A0-x)/4^n/Am;
Y = (A0-y)/4^n/Am;
Z = -X-Y;
E2 = X*Y-Z^2;
E3 = X*Y*Z;
rf = (1-E2/10+E3/14+E2^2/24-3*E2*E3/44)/sqrt(Am);

function [rd]=RD(x,y,z,r)
% RD
% calculates the RD term in the Carlson (1995) method for calculating
% elliptic integrals

% r = 1e-16

if z==0
    error('z value must be nonzero!')
elseif all([x,y]==0)
    error('At most one of the x and y values can be zero!')
end

xm = x;
ym = y;
zm = z;
A0 = (x+y+3*z)/5;
Q = max([abs(A0-x),abs(A0-y),abs(A0-z)])/(r/4)^(1/6);
n = 0;
Am = A0;
S = 0;
while abs(Am)<=Q/(4^n)
    lambdam = sqrt(xm*ym)+sqrt(xm*zm)+sqrt(ym*zm);
    S = S+(1/4^n)/sqrt(zm)/(zm+lambdam);
    Am = (Am+lambdam)/4;
    xm = (xm+lambdam)/4;
    ym = (ym+lambdam)/4;
    zm = (zm+lambdam)/4;
    n = n+1;
end

X = (A0-x)/4^n/Am;
Y = (A0-y)/4^n/Am;
Z = -(X+Y)/3;
E2 = X*Y-6*Z^2;
E3 = (3*X*Y-8*Z^2)*Z;
E4 = 3*(X*Y-Z^2)*Z^2;
E5 = X*Y*Z^3;
rd = (1-3*E2/14+E3/6+9*E2^2/88-3*E4/22-9*E2*E3/52+3*E5/26)/4^n/Am^1.5+3*S;


function [Xs,Ys,Zs,Ws]=AdaptiveGrid(X0,Y0,depth,omegaX,omegaY,omegaZ,...
    ax,ay,az,Cr,Nmax)
% AdaptiveGrid
% Calulates the coordinates and respectives weights of the point CDMs that
% form the solution. This function uses the adaptive algorithm explained in
% the reference paper #1.

% Rotation matrix
Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

% Determine aC, dT, aV, np
aC = max([ax ay az]); % equation 2
[Ztop,Zbot] = EllTopBot(depth,omegaX,omegaY,omegaZ,ax,ay,az);
aV = (Ztop-Zbot)/2; % equation 4
dT = -Ztop;
k1 = dT/(aC*Cr-aV); % equation 6

% Determine aH
normVec = [0 0 1]'; % Unit normal vector of a horizontal plane
normVecR = R'*normVec; % rotated normal vector
Pt = [0 0 0]'; % A point on the plane
PtR = R'*Pt; % A point on the rotated plane
[~,aH] = Ell_Plane_Intersect(ax,ay,az,normVecR(1),normVecR(2),...
    normVecR(3),normVecR'*PtR);

SFa = 0;
SFb = 1;
Nta = 0;
Ntb = Nmax;
itr = 1;
while Nta~=Ntb && abs(SFa-SFb)>1e-3
    
    if itr==1
        SF = 1;
    else
        SF = (SFa+SFb)/2;
    end
    
    np = SF*dT/(2*k1*aV); % equation 8
    nH = SF*dT/(2*k1*aH);
    
    % Determine the partitioning planes: Zpar contains the Z coordinates of
    % the patitioning planes
    Npar0 = floor(log(Zbot/Ztop)/log(1+1/np)); % Initial # of partitions
    if Npar0<=1
        Zpar = [Ztop;Zbot];
    else
        ZparTmp = Ztop*(1+1/np).^(0:Npar0)';
        if ZparTmp(end)>Zbot
            Zpar = [ZparTmp;Zbot];
        else
            Zpar = ZparTmp;
        end
    end
    Npar = numel(Zpar)-1; % Total number of partitioning planes
    Z0 = (Zpar(1:Npar)+Zpar(2:Npar+1))/2; % Planes containing the pCDMs 
    
    % Determine the intersection ellipses
    Ptj = [zeros(2,Npar); Z0(:)'+depth]; % A point on each plane
    PtjR = R'*Ptj; % A point on the rotated plane
    [Pc,ajMax,ajMin,ejMax,ejMin] = Ell_Plane_Intersect(ax,ay,az,...
        normVecR(1),normVecR(2),normVecR(3),normVecR'*PtjR);
    
    sjMax = -Zpar(1:Npar)/nH; % source spacing on the partitioning planes
    njMax = floor(ajMax./sjMax);
    
    sjMin = ajMin./ajMax.*sjMax;
    njMin = ceil(ajMin./sjMin);
    
    Nsp = zeros(Npar,1);
    for kk=1:Npar
        [xs,ys] = meshgrid(-njMax(kk):njMax(kk),-njMin(kk):njMin(kk));
        xs = xs*sjMax(kk);
        ys = ys*sjMin(kk);
        Ind = xs.^2/ajMax(kk)^2+ys.^2/ajMin(kk)^2<1;
        Nsp(kk) = nnz(Ind);
    end
    Nt = sum(Nsp);
    
    if itr==1 && Nt<Nmax
        break
    else
        itr = inf;
    end
    
    if sign(Nt-Nmax)==sign(Nta-Nmax)
        SFa = SF;
        Nta = Nt;
    else
        SFb = SF;
        Ntb = Nt;
    end
end

% Calculate volume fraction of all partitions: Vf
if Ztop==Zbot
    Vf = 1;
else
    Vp = zeros(Npar+1,1); % Volume of the top ellipsoidal caps
    for kk=1:Npar+1
        Vp(kk) = Ell_Cap_Vol(depth,omegaX,omegaY,omegaZ,ax,ay,az,Zpar(kk));
    end
    Vp = diff(Vp); % Volume of the partitions
    Vf = Vp/sum(Vp); % Volume fraction of the partitions
end

Xs = zeros(Nt,1);
Ys = zeros(Nt,1);
Zs = zeros(Nt,1);
Ws = zeros(Nt,1);
Ns0 = 0;
for kk=1:Npar
    
    [xs,ys] = meshgrid(-njMax(kk):njMax(kk),-njMin(kk):njMin(kk));
    xs = xs*sjMax(kk);
    ys = ys*sjMin(kk);
    Ind = xs.^2/ajMax(kk)^2+ys.^2/ajMin(kk)^2<1;
    xs = xs(Ind);
    ys = ys(Ind);
    
    r = [ejMax ejMin normVecR]*[xs(:)';ys(:)';zeros(1,numel(xs))];
    Xr = r(1,:)'+Pc(kk,1);
    Yr = r(2,:)'+Pc(kk,2);
    Zr = r(3,:)'+Pc(kk,3);
    
    r1 = R*[Xr Yr Zr]';
    Xs0 = r1(1,:)'+X0;
    Ys0 = r1(2,:)'+Y0;
    Zs0 = r1(3,:)'-depth;
    
    Nsp = numel(Xs0);
    if Npar>1 && kk==1
        [~,~,~,~,Zs1,Zs2] = Ell_Line_Intersect(X0,Y0,depth,...
            omegaX,omegaY,omegaZ,ax,ay,az,0,0,1,Xs0,Ys0,Zs0);
        Wz = max([Zs1(:) Zs2(:)],[],2)'-Zs0(1);
        Ws0 = Wz(:)+sjMax(kk)/2;
        Ws0 = Vf(kk)*Ws0/sum(Ws0);
    elseif Npar>1 && kk==Npar
        [~,~,~,~,Zs1,Zs2] = Ell_Line_Intersect(X0,Y0,depth,...
            omegaX,omegaY,omegaZ,ax,ay,az,0,0,1,Xs0,Ys0,Zs0);
        Wz = Zs0(1)-min([Zs1(:) Zs2(:)],[],2)';
        Ws0 = Wz(:)+sjMax(kk)/2;
        Ws0 = Vf(kk)*Ws0/sum(Ws0);
    elseif Npar==1
        [~,~,~,~,Zs1,Zs2] = Ell_Line_Intersect(X0,Y0,depth,...
            omegaX,omegaY,omegaZ,ax,ay,az,0,0,1,Xs0,Ys0,Zs0);
        Wz = abs(Zs1-Zs2);
        Ws0 = Wz(:);
        Ws0 = Vf(kk)*Ws0/sum(Ws0);
    else
        Ws0 = Vf(kk)/Nsp;
    end
    
    Xs(Ns0+1:Ns0+Nsp) = Xs0;
    Ys(Ns0+1:Ns0+Nsp) = Ys0;
    Zs(Ns0+1:Ns0+Nsp) = Zs0;
    Ws(Ns0+1:Ns0+Nsp) = Ws0;
    
    Ns0 = Ns0+Nsp;
end


function [Ztop,Zbot,Th,La]=EllTopBot(depth,omegaX,omegaY,omegaZ,ax,ay,az)
% EllTopBot
% Calculates the spherical coordinates and the Z Cartesian coordinates of
% the shallowest (top) and deepest (bottom) points of an ellipsoid. The 
% ellipsoid equation before applying the rotations (omegaX,omegaY,omegaZ) 
% is: x^2/ax^2+y^2/ay^2+z^2/az^2 = 1.

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

La = atan2(R(3,2)*ay,R(3,1)*ax);
Th = atan2(R(3,1)^2*ax^2+R(3,2)^2*ay^2,...
    R(3,3)*az*sqrt(R(3,1)^2*ax^2+R(3,2)^2*ay^2));

if La<0
    La = La+2*pi;
end

if Th<0
    Th = Th+pi;
elseif Th>pi
    Th = Th-pi;
end

% Radian to degree
La = La*180/pi;
Th = Th*180/pi;

Z = R(3,1)*ax*sind(Th)*cosd(La)+...
    R(3,2)*ay*sind(Th)*sind(La)+R(3,3)*az*cosd(Th);
Z = [Z -Z]-depth;
Ztop = max(Z);
Zbot = min(Z);


function [Vtop]=Ell_Cap_Vol(depth,omegaX,omegaY,omegaZ,ax,ay,az,Z0)
% Ell_Cap_Vol
% calculates the ellipsoidal cap volume formed by the intersection of the
% horizontal plane z = Z0 and an ellipsoid. The ellipsoid equation before 
% applying the rotations (omegaX,omegaY,omegaZ) is: 
% x^2/ax^2+y^2/ay^2+z^2/az^2 = 1 and the plane equation is: Ax+By+Cz = D.


Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

normVec = [0 0 1]'; % Unit normal vector of the plane
normVecR = R'*normVec; % rotated normal vector
Pt = [0 0 Z0+depth]'; % A point on the plane
PtR = R'*Pt; % A point on the rotated plane

A = normVecR(1);
B = normVecR(2);
C = normVecR(3);
D = normVecR'*PtR;
dn = D/sqrt((A*ax)^2+(B*ay)^2+(C*az)^2);
Vtop = pi*ax*ay*az/3*(1-dn)^2*(2+dn);


function [Pc,aMaj,aMin,eMaj,eMin]=Ell_Plane_Intersect(ax,ay,az,a,b,c,d)
% Ell_Plane_Intersect
% Finds the ellipse formed by the intersection of an ellipsoid and a plane.
% The ellipsoid equation before applying the rotations 
% (omegaX,omegaY,omegaZ) is: x^2/ax^2+y^2/ay^2+z^2/az^2 = 1 and the plane 
% equation is: ax+by+cz = d. Here only "d" can be vector.
% 
% This function is based on Reference #6: Klein (2012). 

% Find a vector that is normal to [a b c]
[mx,Ind] = max(abs([a b c]));
if Ind==1
    r = cross([a b c],[0 mx mx])';
elseif Ind==2
    r = cross([a b c],[mx 0 mx])';
elseif Ind==3
    r = cross([a b c],[mx mx 0])';
end
r = r/sqrt(r'*r);
normVec = [a b c]'/sqrt(a^2+b^2+c^2);
s = cross(normVec,r);
D1 = diag(1./[ax ay az]);
w = atan2(2*(D1*r)'*(D1*s),((D1*r)'*(D1*r)-(D1*s)'*(D1*s)))/2;
e1 = cos(w)*r+sin(w)*s;
e2 = -sin(w)*r+cos(w)*s;

d = d(:);
kappa = d/sqrt(a^2+b^2+c^2); % distance from the origin
dn = kappa.^2*(a^2+b^2+c^2)/((a*ax)^2+(b*ay)^2+(c*az)^2);

beta1 = (D1*e1)'*(D1*e1);
beta2 = (D1*e2)'*(D1*e2);
a1 = sqrt((1-dn)/beta1);
a2 = sqrt((1-dn)/beta2);
if a1(1)>a2(1)
    aMaj = a1;
    eMaj = e1;
    aMin = a2;
    eMin = e2;
else
    aMaj = a2;
    eMaj = e2;
    aMin = a1;
    eMin = e1;
end

Ck = kappa*(a^2+b^2+c^2)/((a*ax)^2+(b*ay)^2+(c*az)^2)/sqrt(a^2+b^2+c^2);
Pc = [Ck*ax^2*a, Ck*ay^2*b, Ck*az^2*c];


function [Xs1,Xs2,Ys1,Ys2,Zs1,Zs2]=Ell_Line_Intersect(X0,Y0,depth,...
    omegaX,omegaY,omegaZ,ax,ay,az,a0,b0,c0,x0,y0,z0)
% Ell_Line_Intersect
% Finds the intersection points of a line and an ellipsoid. The code can 
% handle one ellipsoid and multiple lines. The ellipsoid equation before 
% the rotations are applied is: x^2/ax^2+y^2/ay^2+z^2/az^2 = 1. The
% genaral equation of the lines are: x = ka0+x0 ,y = kb0+y0, z = kc0+z0.

Rx = [1 0 0;0 cosd(omegaX) sind(omegaX);0 -sind(omegaX) cosd(omegaX)];
Ry = [cosd(omegaY) 0 -sind(omegaY);0 1 0;sind(omegaY) 0 cosd(omegaY)];
Rz = [cosd(omegaZ) sind(omegaZ) 0;-sind(omegaZ) cosd(omegaZ) 0;0 0 1];
R = Rz*Ry*Rx;

if numel(a0)==1 && numel(x0)>1
    a0 = repmat(a0,size(x0));
    b0 = repmat(b0,size(x0));
    c0 = repmat(c0,size(x0));
end

unitVec = [a0(:) b0(:) c0(:)]'; % Unit vector along the line
P0 = [x0(:)-X0 y0(:)-Y0 z0(:)+depth]'; % A point on the line

unitVecR = R'*unitVec; % rotated unit vector
P0R = R'*P0; % rotated point
x1 = P0R(1,:);
y1 = P0R(2,:);
z1 = P0R(3,:);

a = unitVecR(1,:);
b = unitVecR(2,:);
c = unitVecR(3,:);

% Intersection of a generic line and a standard ellipsoid
A = a.^2*ay^2*az.^2 + b.^2*ax^2*az^2 + c.^2*ax^2*ay^2;
B = 2*(a.*x1*ay^2*az^2 + b.*y1*ax^2*az^2 + c.*z1*ax^2*ay^2);
C = x1.^2*ay^2*az^2 + y1.^2*ax^2*az^2 + z1.^2*ax^2*ay^2 - ax^2*ay^2*az^2;
[k1,k2] = SolveQuadraticVec(A,B,C);

xs1 =  k1.*a+x1;
xs2 =  k2.*a+x1;
ys1 =  k1.*b+y1;
ys2 =  k2.*b+y1;
zs1 =  k1.*c+z1;
zs2 =  k2.*c+z1;

r1 = R*[xs1(:) ys1(:) zs1(:)]';
r2 = R*[xs2(:) ys2(:) zs2(:)]';
Xs1 = r1(1,:)'+X0;
Ys1 = r1(2,:)'+Y0;
Zs1 = r1(3,:)'-depth;
Xs2 = r2(1,:)'+X0;
Ys2 = r2(2,:)'+Y0;
Zs2 = r2(3,:)'-depth;


function [x1,x2]=SolveQuadraticVec(A,B,C)
% SolveQuadraticVec
% Finds the roots of A(k)*x^2+B(k)*x+C(k) = 0, where k=1,...,N and N is the
% total number of elements of A. The matrices A, B and C must have the same
% size.

D = B.^2-4*A.*C;

x1 = zeros(size(A));
x2 = x1;

Ind = B>=0;
x1(Ind) = (-B(Ind)-sqrt(D(Ind)))/2./A(Ind);
x2(Ind) = 2*C(Ind)./(-B(Ind)-sqrt(D(Ind)));
x1(~Ind) = 2*C(~Ind)./(-B(~Ind)+sqrt(D(~Ind)));
x2(~Ind) = (-B(~Ind)+sqrt(D(~Ind)))/2./A(~Ind);

% In this particular problem D cannot be negative!
x1(D<0) = 0;
x2(D<0) = 0;