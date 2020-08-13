% Numerical integral route to the Redfield relaxation superopera-
% tor. Syntax:
%
%                R=ngce(spin_system,H0,H1,traj_pts,dt)
%
% Parameters:
%
%  H0 - static laboratory frame Hamiltonian commutation su-
%       peroperator acting in the background, a matrix
%
%  H1 - stochastic part (zero mean) of the laboratory frame
%       Hamiltonian commutation superoperator, a cell array
%       of matrices for each point in the MD trajectory.
%
%  traj_pts - number of points in the trajectory segment
%
%  dt - time step of the MD trajectory, seconds
%           
% Outputs:
%
%  R  - laboratory frame relaxation superoperator
%
% Notes: enough trajectory points (number of segments x traj_pts) must be 
% present to converge the ensemble averages and Redfield's integral.
%
% i.kuprov@soton.ac.uk
%
% <http://spindynamics.org/wiki/index.php?title=ngce.m>
% modified by J. Prestegard 11/10/19 to run with coordinate sampling of
% methyl protons plus one - removed sampling checks lines 67-77 for testing
% purposes - also chnaged norm(H0,2) to norm(H0,1) 02/03/20
% Changed to return R for one segment only 02/24/20 - averaging done in ngce_test 
% Removed error checking for total trajectory length

function R_gce=ngce_IK_JHP_3(spin_system,H0,H1,traj_pts,dt)

% Check consistency
grumble(spin_system,H0,H1,dt);

% Coherent dynamics timescale
timescale=2*pi/norm(H0,1);

% Number of points under tau integral
n_tau_int_steps=traj_pts;


% Trajectory duration
traj_dur=dt*(traj_pts-1);

% Print timing diagnostics
report(spin_system,' ');
report(spin_system,['H0 period:                       ' num2str(timescale,'%12.5e') ' seconds']);
report(spin_system,['Trajectory step length:          ' num2str(dt,'%12.5e') ' seconds']);
report(spin_system,['Total trajectory length:         ' num2str(traj_dur,'%12.5e') ' seconds']);
report(spin_system,['Total trajectory points:         ' num2str(traj_pts)]);
%report(spin_system,['User estimate for tau_c:         ' num2str(tau_est,'%12.5e') ' seconds']);

% Enforce sufficient H0 sampling
report(spin_system,['Trajectory points per H0 period: ' num2str(timescale/dt)]);
if timescale/dt<50
    error('insufficient H0 dynamics sampling, reduce trajectory time step.');
end

% Enforce sufficient tau_c sampling
%report(spin_system,['Trajectory points per tau_est:   ' num2str(tau_est/dt)]);
if traj_pts<25
    error('insufficient correlation function sampling, reduce trajectory time step.');
end

% Enforce sufficient trajectory duration
%report(spin_system,['Trajectory duration / tau_est:   ' num2str(dt*(traj_pts-1)/tau_est)]);
%if traj_pts<200
    %error('insufficient ensemble average sampling, increase trajectory duration.');
%end

% Get propagators for tau integrals
report(spin_system,'computing H0 exponentials...');
P=cell(n_tau_int_steps,1);

P_dt=expm(-1i*H0*dt); 
P{1}=eye(size(H0));
for n=2:n_tau_int_steps
    P{n}=P_dt*P{n-1};
end

% Preallocate the answer
R_gce=zeros(size(H0));

% Compute trajectory integrals
report(spin_system,'computing Redfield''s integral...'); 
% Preallocate current integral
F=zeros(size(H0));
f_curr=zeros(size(H0));
    
% get H1 trajecory segment for start point
    % Trapezium rule for tau integral
    f_curr=P{1}*H1{1}*P{1}';
    for tau=1:traj_pts-1
        f_next=P{1+tau}*H1{1+tau}*P{1+tau}';
        F=F+(f_curr+f_next)/2;
        f_curr=f_next;
    end
   
R_gce=-dt*H1{1}*F;
clear('H1');
end

% Consistency enforcement
function grumble(spin_system,H0,H1,dt)
if ~isnumeric(H0)
    error('H0 must be a matrix.');
end
if ~iscell(H1)
    error('H1 must be a cell array of matrices.');
end
if (~isnumeric(dt))||(~isreal(dt))||(dt<=0)
    error('dt must be a positive real number.');
end
%if (~isnumeric(tau_est))||(~isreal(tau_est))||(tau_est<=0)
    %error('tau_est must be a positive real number.');
%end
end

% "To engage in hand-to-hand combat, a trooper must have lost
% somewhere on the battlefield his rifle, pistol, knife, belt, 
% shovel, vest, and helmet. He must also have found a comple-
% tely clean spot with not a single stick or stone on it. And
% only then can he engage into ferocious hand-to-hand combat
% with a similar dumbass."
%
% A Russian hand-to-hand combat instructor

