function PWP(met_input_file, profile_input_file, pwp_output_file)

%  function PWP(met_input_file, profile_input_file, pwp_output_file)

%  MATLAB version 1.4 (8 Oct 2001)

%  This main program is used to drive the pwp subroutine (pwpgo)
%  with a daily cycle of heating, or, an arbitrary (observed) time
%  series of surface flux data.

%  Written by Jim Price, april 27, 1989. Please direct any
%  comments or questions to Jim Price, WHOI, Woods Hole, MA
%  02543, USA, tel. 508-289-2526.

%  Version 2:  changed the gradient Ri relaxation method.
%  Version 3:  clarified the bulk Ri relaxation.
%  Version 4:  fixed an error that disabled background diffusion.

%  Last editted on 20 Jan, 1999 by JFP.

%  Translated for MATLAB 5.X by Peter Lazarevich and Scott Stoermer,
%  GSO/URI, 11 July 2001.


%  -- Option to display output during run. --

% BH - IGNORE THIS
diagnostics = 0;

%  -- Initialize the user-defined variables. --

% BH - IGNORE THIS
global dt dz days depth dt_save lat g cpw rb rg rkz beta1 beta2 udrag

%  dt,			time-step increment (seconds)
%  dz,			depth increment (meters)
%  days,  	the number of days to run
%  depth,		the depth to run
%  dt_save,	time-step increment for saving to file (multiples of dt)
%  lat,  		latitude (degrees)
%  g,				gravity (9.8 m/s^2)
%  cpw,			specific heat of water (4183.3 J/kgC)
%  rb,			critical bulk richardson number (0.65)
%  rg,			critical gradient richardson number (0.25)
%  rkz,			background diffusion (0)
%  beta1, 	longwave extinction coefficient (0.6 m)
%  beta2, 	shortwave extinction coefficient (20 m)

% TRANSLATE USING <- INSTEAD OF =. REMOVE ; AT END OF LINES (WEIRD MATLAB THING)
dt			= 900;
dz			= 10;
days 		= 300
depth		= 1000
dt_save	= 4;
lat 		= 60
g				= 9.8;
cpw			= 4183.3;
rb			= 0.65;
rg			= 0.25;
rkz			= 0;
beta1 	= 0.6;
beta2 	= 20;

%  -- Initialize additional variables. ---

% IGNORE
global f ucon

%  f,			Coriolis parameter
%  ucon,	coefficient of inertial-internal wave dissipation

% TRANSLATE AS IS (WITH CHANGING = TO <-)
f = 2*7.29E-5*sin(lat*pi/180);

%  -- Load the air/sea flux data. --

% IGNORE
disp(['loading ' met_input_file])

% ASSUME THAT WE ARE GOING TO BE READING IN OUR MET DATA FROM A CSV FILE. USE read_csv FROM TIDYVERSE AND ASSIGN TO An object "met"
load(met_input_file);

%  -- Interpolate the air/sea flux variables at dt resolution. --

%  Set nmet equal to the number of time increments using a resolution of dt.

% CHANGE . TO $ FOR REFERENCING COLUMN OF met. CHANGE () TO [] FOR REFERENCING FIRST VALUE OF TIME
nmet 	= days*8.64E4/dt;
time	= met.time(1)+(0:(nmet-1))*dt/8.64E4;

% Check record length
% IGNORE FOR NOW
if time(end) > met.time(end)
	time = time(time<met.time(end));
	nmet = length(time);
	disp(['Met input shorter than # of days selected, truncating run to ' num2str(nmet*dt/8.64E4) ' day(s)'])
	pause
end

% Check the time-resolution of the inertial period
% IGNORE FOR NOW
if dt > 1/10*2*pi/f
	ans = input('Time step, dt, too large to accurately resolve the inertial period. Is this okay? (y/n)','s');
	if ans == 'n'
		disp(['Please restart PWP.m with a new dt <= ' num2str(1/(10*f))])
		break
	end
end

% REPLACE interp1 WITH approx() - USE THE HELP FILE TO WORK OUT HOW THAT WORKS AS WELL AS USING THE NEAREST POINT INTERPOLATION LIKE IS USED HERE
% AGAIN, CHANGE . TO $ FOR SELECTING A COLUMN IN THE DATA. NOT SURE WHY THEY USE floor HERE. I'D REMOVE.
qi			= interp1(met.time,met.sw,floor(time),'nearest');
qo			= interp1(met.time,met.lw+met.qlat+met.qsens,floor(time),'nearest');
tx			= interp1(met.time,met.tx,floor(time),'nearest');
ty			= interp1(met.time,met.ty,floor(time),'nearest');
precip	= interp1(met.time,met.precip,floor(time),'nearest');

% IGNORE
disp('loading complete')

%  -- Load initial t,s profile data. --

% IGNORE
disp(['loading ' profile_input_file])

% AGAIN ASSUME WE ARE GOING TO READ IN THE
load(profile_input_file);

%  -- Interpolate the profile variables at dz resolution. --

% IGNORE
global nz z t s d

%  Set nz equal to the number of depth increments + 1 using a resolution of dz.
% SHOULD BE SIMPLE TRANSLATION
nz	= 1+depth/dz;
z		= ((0:nz-1)*dz)';

% Check record length
% IGNORE FOR NOW
if z(end) > profile.z(end)
	z = z(z<=profile.z(end));
	nz = length(z);
	disp(['Profile input shorter than depth selected, truncating to ' num2str(z(end)) ' meters'])
	pause
end

% Check the depth-resolution of the profile file
% IGNORE FOR NOW
profile_increment = (profile.z(end)-profile.z(1))/(length(profile.z)-1);
if dz < profile_increment/5
	ans = input('Depth increment, dz, is much smaller than profile resolution. Is this okay? (y/n)','s');
	if ans == 'n'
		disp(['Please restart PWP.m with a new dz >= ' num2str(profile_increment/5)])
		break
	end
end

% AS ABOVE, REPLACE WITH approx()
t	= interp1(profile.z,profile.t,z);
s	= interp1(profile.z,profile.s,z);
d	= sw_dens0(s,t);

% IGNORE
disp('loading complete')

%  -- Interpolate evaporation minus precipitaion at dt resolution. --

% AS ABOVE
evap 	= (0.03456/(86400*1000))*interp1(met.time,met.qlat,floor(time),'nearest');
emp		= evap - precip;


% HAVE THIS AS THE FIRST GOAL FOR TRANSLATION - NEXT STEP IS TO GO TO pwpgo FUNCTION. LET ME KNOW WHEN YOU WANT COMMENTS IN THAT.


%  -- Initialize additional profile variables at dz resolution. --


global u v absrb

%  u and v, 	east and north current
%  absrb,			absorbtion fraction

% NEED TO CHANGE zeros TO rep(0,nz)
% NEED TO TRANSLATE absorb FUNCTION. SEE BELOW.
u			= zeros(nz,1);
v			= zeros(nz,1);
absrb = absorb(beta1,beta2);

%  -- Specify a simple "background" diffusion to be applied to the profiles. --

dstab = dt*rkz/dz^2;

if dstab > 0.5
	disp('Warning, this value of rkz will be unstable')
end

%  -- Define the variables to be saved. --

pwp_output.dt 			= dt;
pwp_output.dz 			= dz;
pwp_output.lat 			= lat;
pwp_output.z				= z;
pwp_output.time			= [];
pwp_output.t				= [];
pwp_output.s				= [];
pwp_output.d				= [];
pwp_output.u				= [];
pwp_output.v				= [];

%  -- Step through the PWP model. --

disp(['STATUS (out of ' int2str(nmet) ' steps):'])

for m = 1:nmet

% LETS START WORKING ON THIS FUNCTION - SEE BELOW
	pwpgo(qi(m),qo(m),emp(m),tx(m),ty(m),m);

%  Apply a "background" diffusion if rkz is non-zero.

	if rkz > 0
		diffus(dstab,t);
		diffus(dstab,s);
		d = sw_dens0(s,t);
		diffus(dstab,u);
		diffus(dstab,v);
	end

%  -- Store variables. --

	if mod(m-1,dt_save) == 0
		pwp_output.time(:,end+1)	= time(m);
		pwp_output.t(:,end+1)			= t;
		pwp_output.s(:,end+1)			= s;
    pwp_output.d(:,end+1)			= d;
    pwp_output.u(:,end+1)			= u;
    pwp_output.v(:,end+1)			= v;
	end

%  -- Save variables. --

	if mod(m, 100) == 0
		disp([int2str(m), ' (' sprintf('%2.1f',100*m/nmet), '%)'])
		eval(['save ' pwp_output_file ' pwp_output'])
	end

%  -- Diagnostics. --

	if ~exist('x')
		x = u(1)*dt;
		y = v(1)*dt;
		U = mean(u)*depth;
		V = mean(v)*depth;
	else
		x = [x x(end)+u(1)*dt];
		y = [y y(end)+v(1)*dt];
		U = [U mean(u)*depth];
		V = [V mean(v)*depth];
	end

	if diagnostics

		figure(1)

		subplot(1,3,1)
		plot(d,z,'k')
		set(gca,'ydir','reverse')
		hold on
		plot(d,z,'k.')
		hold off
		grid on
		title('Density')
		ylabel('Depth (m)')

		subplot(1,3,2)
		plot(t,z,'k')
		set(gca,'ydir','reverse')
		hold on
		plot(t,z,'k.')
		hold off
		grid on
		title('Temp.')

		subplot(1,3,3)
		plot(sqrt(u.^2+v.^2),z,'k')
		set(gca,'ydir','reverse')
		hold on
		plot(sqrt(u.^2+v.^2),z,'k.')
		hold off
		grid on
		title('Speed (m/s)')

		figure(2)

		subplot(1,2,1)
		plot(x,y,'k')
		hold on
		plot(x,y,'k.')
		plot(x(end),y(end),'ko')
		hold off
		grid on
		title('Trajectroy at surface')
		ylabel('y (m)')
		xlabel('x (m)')

		subplot(1,2,2)
		compass(U(end),V(end),'k')
		hold on
		plot(U,V,'k.')
		hold off
		grid on
		title('Net transport (m^3/s)')

		figure(3)

		subplot(2,1,1)
		plot(time(1:m),qi(1:m),'k*')
		hold on
		plot(time(1:m),qo(1:m),'ko')
		hold off
		grid on
		title('qi(*) and qo(o)')
		ylabel('Heat flux (W/m^2)')

		subplot(2,1,2)
		plot(time(1:m),tx(1:m),'k*')
		hold on
		plot(time(1:m),ty(1:m),'ko')
		hold off
		grid on
		title('tx(*) and ty(o)')
		ylabel('Wind stress (N/m^2)')
		xlabel('Time (days)')
		pause(0.2)

	end

end

eval(['save ' pwp_output_file ' pwp_output'])
disp(['Model Run Completed - Selected variables saved in file: ' pwp_output_file])

pwp_output

% -------------------------------------------------------------------------

% CREATING A FUNCTION - SEE CHAPTER 19 OF R4DS (https://r4ds.had.co.nz/functions.html) FOR SYNTAX
function pwpgo(qi,qo,emp,tx,ty,m)

%  This subroutine is an implementation of the Price, Weller,
%  Pinkel upper ocean model described in JGR 91, C7 8411-8427
%  July 15, 1986 (PWP). This version was coded and documented by
%  Jim Price in April 1989.

%  Edited on 20 September 1993 by JFP to allow for a critical
%  gradient Richardson number other than 1/4, and to implement a
%  different and a priori better means of achieving convergence
%  of gradient ri mixing (see subroutine stir for details).
%  The major difference is that the revised scheme gives a more
%  smoothly varying mixed layer depth over a diurnal cycle.
%  Edited on 14 December 1998 by JFP to clarify the bulk Ri
%  relaxation.

%  This model also implements an energy budget form of an
%  entrainment parameterization that is very similar to that
%  described in Price, Mooers, Van Leer, JPO 8, 4, 582-599 (and
%  references therein). This part of the model should be
%  treated as developmental only, as there are several features
%  that are arbitrary to this model (i.e., the depth of the ml
%  during times of heating can be the grid interval, dz). To
%  use this parameterization set the bulk Richardson number to
%  zero, and set em1, or em2, or em3 to non-zero. The gradient
%  Richardson number can be zero or not.

% MATLAB DOES THIS THING WHERE YOU CAN SET GLOBAL VARIABLES - I.E. VARIABLES THAT THE FUNCTION CAN SEE OUTSIDE OF THE FUNCTION.
% THIS IS NOT USUALLY THE CASE FOR MOST LANGUAGES AS IT CAN LEAD TO ISSUES. WE WILL NEED TO FIND A WAY TO HAVE ALL THESE THINGS BE INPUTS INTO THE FUNCTION WE WRITE.
% IF THE FUNCTION SIDE OF THINGS IS CONFUSING, FOCUS ON THE TRANSLATION FIRST AND WE CAN TALK ABOUT THE FUNCTIONALIZATION LATER
global dt dz days depth dt_save lat g cpw rb rg beta1 beta2 udrag
global nz z t s d
global u v absrb
global f ucon

%  Apply heat and fresh water fluxes to the top most grid cell.

% CHANGE INDEXING TO BE [] RATHER THAN (). CHANGE ./ FOR /
t(1) = t(1)+(qi*absrb(1)-qo)*dt./(dz*d(1)*cpw);
s(1) = s(1)/(1-emp*dt/dz);

%  Absorb solar radiation at depth.
% AGAIN, INDEXING AND ./
t(2:nz) = t(2:nz)+qi*absrb(2:nz)*dt./(dz*d(2:nz)*cpw);

%  Compute the density, and relieve static instability, if it occurs.
% USE swSigma AS ABOVE
d = sw_dens0(s,t);

% NEED TO WRITE THIS FUNCTION, LEAVE FOR NOW OF HAVE A LOOK AT IT.
remove_si;

%  At this point the density proifile should be statically stable.

%  Find the index of the surfacd mixed-layer right after the heat/salt fluxes.

% FIND EQUVELENT IN R IS which(). min AND diff are the same.
ml_index = min(find(diff(d)>1E-4));

% AGAIN A CHECK - IGNORE FOR NOW
if isempty(ml_index)
	error_text = ['Error reached in PWP.m (line 346): the mixed layer is too deep!' sprintf('\n') ...
	'Possible reasons:' sprintf('\n') ...
	'1) tempearature inversion at surface is too strong' sprintf('\n') ...
	'2) wind is too strong' sprintf('\n') ...
	'3) evaporation is too great'];
	error([ sprintf('\n') error_text])
end

%  Get the depth of the surfacd mixed-layer.
% INDEX USING []
ml_depth = z(ml_index+1);

%  Time step the momentum equation.

%  Rotate the current throughout the water column through an
%  angle equal to inertial rotation for half a time step.

ang = -f*dt/2;

% CAN REWRITE THIS FUNCTION BELOW
rot(ang);

%  Apply the wind stress to the mixed layer as it now exists.

% INDEXING BUT ELSE ALL GOOD
du = (tx/(ml_depth*d(1)))*dt;
dv = (ty/(ml_depth*d(1)))*dt;
u(1:ml_index) = u(1:ml_index)+du;
v(1:ml_index) = v(1:ml_index)+dv;

%  Apply drag to the current (this is a horrible parameterization of
%  inertial-internal wave dispersion).
% DONT REALLY KNOW WHAT THIS DOES, BUT CAN CHANGE TO IF STATEMENT IN R
if ucon > 1E-10
	u = u*(1-dt*ucon);
	v = v*(1-dt*ucon);
end

%  Rotate another half time step.
% AS ABOVE
rot(ang);

%  Finished with the momentum equation for this time step.

%  Do the bulk Richardson number instability form of mixing (as in PWP).

% WILL LOOK INTO THIS FUNCTION AND THE ONE BELOW AND THINK ABOUT BEST APPROACH.
% HONESTLY NONE OF THEM ARE LOOKING LIKE THERES TOO MUCH TO TRANSLATE IN TERMS OF FINDING EQUIVELENT FUNCTIONS
% MOSTLY JUST CHANGING HOW THINGS ARE INDEXED.
if rb > 1E-5
	bulk_mix(ml_index)
end

%  Do the gradient Richardson number instability form of mixing.

if rg > 0
	grad_mix;
end


% -------------------------------------------------------------------------

% TEST COMMENT
function bulk_mix(ml_index)

global g rb
global nz z d
global u v

rvc = rb;
for j = ml_index+1:nz
	h 	= z(j);
	dd 	= (d(j)-d(1))/d(1);
	dv 	= (u(j)-u(1))^2+(v(j)-v(1))^2;
	if dv == 0
		rv = Inf;
	else
		rv = g*h*dd/dv;
	end
	if rv > rvc
		break
	else
		mix5(j);
	end
end


% -------------------------------------------------------------------------


function grad_mix

%  This function performs the gradeint Richardson Number relaxation
%  by mixing adjacent cells just enough to bring them to a new
%  Richardson Number.

global dz g rg
global nz z t s d
global u v

rc 	= rg;

%  Compute the gradeint Richardson Number, taking care to avoid dividing by
%  zero in the mixed layer.  The numerical values of the minimum allowable
%  density and velocity differences are entirely arbitrary, and should not
%  effect the calculations (except that on some occasions they evidnetly have!)

j1 = 1;
j2 = nz-1;

while 1
	for j = j1:j2
		if j <= 0
			keyboard
		end
		dd = (d(j+1)-d(j))/d(j);
		dv = (u(j+1)-u(j))^2+(v(j+1)-v(j))^2;
		if dv == 0
			r(j) = Inf;
		else
			r(j) = g*dz*dd/dv;
		end
	end

	%  Find the smallest value of r in profile

	rs = min(r);
	js = min(find(r==rs));

	%  Check to see whether the smallest r is critical or not.

	if rs > rc
		return
	end

	%  Mix the cells js and js+1 that had the smallest Richardson Number

	stir(rc,rs,js);

	%  Recompute the Richardson Number over the part of the profile that has changed

	j1 = js-2;
	if j1 < 1
		 j1 = 1;
	end
	j2 = js+2;
	if j2 > nz-1
		 j2 = nz-1;
	end
end


% -------------------------------------------------------------------------


function a = stir(rc,r,j)

%  This subroutine mixes cells j and j+1 just enough so that
%  the Richardson number after the mixing is brought up to
%  the value rnew. In order to have this mixing process
%  converge, rnew must exceed the critical value of the
%  richardson number where mixing is presumed to start. If
%  r critical = rc = 0.25 (the nominal value), and r = 0.20, then
%  rnew = 0.3 would be reasonable. If r were smaller, then a
%  larger value of rnew - rc is used to hasten convergence.

%  This subroutine was modified by JFP in Sep 93 to allow for an
%  aribtrary rc and to achieve faster convergence.

global t s d u v

rcon 			= 0.02+(rc-r)/2;
rnew 			= rc+rcon/5;
f 				= 1-r/rnew;
dt				= (t(j+1)-t(j))*f/2;
t(j+1)		= t(j+1)-dt;
t(j) 			= t(j)+dt;
ds				= (s(j+1)-s(j))*f/2;
s(j+1)		= s(j+1)-ds;
s(j) 			= s(j)+ds;
d(j:j+1)	= sw_dens0(s(j:j+1),t(j:j+1));
du				= (u(j+1)-u(j))*f/2;
u(j+1)		= u(j+1)-du;
u(j) 			= u(j)+du;
dv				= (v(j+1)-v(j))*f/2;
v(j+1)		= v(j+1)-dv;
v(j) 			= v(j)+dv;


% -------------------------------------------------------------------------


function mix5(j)

%  This subroutine mixes the arrays t, s, u, v down to level j.

global t s d u v

t(1:j) = mean(t(1:j));
s(1:j) = mean(s(1:j));
d(1:j) = sw_dens0(s(1:j),t(1:j));
u(1:j) = mean(u(1:j));
v(1:j) = mean(v(1:j));


% -------------------------------------------------------------------------


function rot(ang)

%  This subroutine rotates the vector (u,v) through an angle, ang

global u v

r = (u+i*v)*exp(i*ang);
u = real(r);
v = imag(r);


% -------------------------------------------------------------------------


function remove_si

%  Find and relieve static instability that may occur in the
%  density array d. This simulates free convection.
%  ml_index is the index of the depth of the surface mixed layer after adjustment,

global d

while 1
	ml_index = min(find(diff(d)<0));
	if isempty(ml_index)
		break
	end
	mix5(ml_index+1);
end


% -------------------------------------------------------------------------


function absrb = absorb(beta1,beta2)

%  Compute solar radiation absorption profile. This
%  subroutine assumes two wavelengths, and a double
%  exponential depth dependence for absorption.

%  Subscript 1 is for red, non-penetrating light, and
%  2 is for blue, penetrating light. rs1 is the fraction
%  assumed to be red.

global nz dz

rs1 = 0.6;
rs2 = 1.0-rs1;
absrb = zeros(nz,1);
z1 = (0:nz-1)*dz;
z2 = z1 + dz;
z1b1 = z1/beta1;
z2b1 = z2/beta1;
z1b2 = z1/beta2;
z2b2 = z2/beta2;
absrb = (rs1*(exp(-z1b1)-exp(-z2b1))+rs2*(exp(-z1b2)-exp(-z2b2)))';

% -------------------------------------------------------------------------


function a = diffus(dstab,a)

%  This subroutine applies a simple diffusion
%  operation to the array a. It leaves the endpoints
%  unchanged (assumes nothing about the
%  boundary conditions).

global nz

a(2:nz-1) = a(2:nz-1)+dstab*(a(1:nz-2)-2*a(2:nz-1)+a(3:nz));


% -------------------------------------------------------------------------


function dens = sw_dens0(S,T)

% SW_DENS0   Denisty of sea water at atmospheric pressure
%=========================================================================
% SW_DENS0  $Revision: 1.3 $  $Date: 1994/10/10 04:54:09 $
%           Copyright (C) CSIRO, Phil Morgan 1992
%
% USAGE:  dens0 = sw_dens0(S,T)
%
% DESCRIPTION:
%    Density of Sea Water at atmospheric pressure using
%    UNESCO 1983 (EOS 1980) polynomial.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   dens0 = density  [kg/m^3] of salt water with properties S,T,
%           P=0 (0 db gauge pressure)
%
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%     Unesco 1983. Algorithms for computation of fundamental properties of
%     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%
%     Millero, F.J. and  Poisson, A.
%     International one-atmosphere equation of state of seawater.
%     Deep-Sea Res. 1981. Vol28A(6) pp625-629.
%=========================================================================

% CALLER: general purpose, sw_dens.m
% CALLEE: sw_smow.m

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=2
   error('sw_dens0.m: Must pass 2 parameters')
end %if

[mS,nS] = size(S);
[mT,nT] = size(T);

if (mS~=mT) | (nS~=nT)
   error('sw_dens0.m: S,T inputs must have the same dimensions')
end %if

Transpose = 0;
if mS == 1  % a row vector
  S = S(:);
  T = T(:);
  Transpose = 1;
end %if

%----------------------
% DEFINE CONSTANTS
%----------------------
%     UNESCO 1983 eqn(13) p17.

b0 =  8.24493e-1;
b1 = -4.0899e-3;
b2 =  7.6438e-5;
b3 = -8.2467e-7;
b4 =  5.3875e-9;

c0 = -5.72466e-3;
c1 = +1.0227e-4;
c2 = -1.6546e-6;

d0 = 4.8314e-4;

%$$$ dens = sw_smow(T) + (b0 + b1*T + b2*T.^2 + b3*T.^3 + b4*T.^4).*S  ...
%$$$                    + (c0 + c1*T + c2*T.^2).*S.*sqrt(S) + d0*S.^2;

dens = sw_smow(T) + (b0 + (b1 + (b2 + (b3 + b4*T).*T).*T).*T).*S  ...
                   + (c0 + (c1 + c2*T).*T).*S.*sqrt(S) + d0*S.^2;

if Transpose
  dens = dens';
end %if

return


% -------------------------------------------------------------------------


function dens = sw_smow(T)

% SW_SMOW    Denisty of standard mean ocean water (pure water)
%=========================================================================
% SW_SMOW  $Revision: 1.3 $  $Date: 1994/10/10 05:51:46 $
%          Copyright (C) CSIRO, Phil Morgan 1992.
%
% USAGE:  dens = sw_smow(T)
%
% DESCRIPTION:
%    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
%
% INPUT:
%   T = temperature [degree C (IPTS-68)]
%
% OUTPUT:
%   dens = density  [kg/m^3]
%
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%     Unesco 1983. Algorithms for computation of fundamental properties of
%     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%     UNESCO 1983 p17  Eqn(14)
%
%     Millero, F.J & Poisson, A.
%     INternational one-atmosphere equation of state for seawater.
%     Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
%=========================================================================

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
% TEST INPUTS
if nargin ~= 1
   error('sw_smow.m: Only one input argument allowed')
end %if

Transpose = 0;
[mT,nT] = size(T);
if mT == 1 % a row vector
   T = T(:);
   Tranpose = 1;
end %if

%----------------------
% DEFINE CONSTANTS
%----------------------
a0 = 999.842594;
a1 =   6.793952e-2;
a2 =  -9.095290e-3;
a3 =   1.001685e-4;
a4 =  -1.120083e-6;
a5 =   6.536332e-9;

dens = a0 + (a1 + (a2 + (a3 + (a4 + a5*T).*T).*T).*T).*T;

if Transpose
  dens = dens';
end %if

return
