% NOE build-up curves using a stored relaxation matrix.  Initial model: 
% three methyl protons and a remote proton. 
% Based on a 2_spin_NOE script by I. Kuprove i.kuprov@soton.ac.uk
% modified by J Prestegard 07/09/20
% modified by J Prestegard 08/03/20 to correct for numerical integration
% errors as suggested by I Kuprov
%
% function noe_build_up_from_stored_R_matix()
%
%  reg - optional overall relaxation rate, this is added to 
%        every eigenvalue of the resulting matrix to prevent
%        very small relaxation rates (e.g. singlets) from 
%        jumping into positive due to integration accuracy
%        limits and then causing problems
%   
reg = 0.2;
%
% function noe_build_up_from_stored_R_matix()
%
% Set the spin system
sys.isotopes={'1H','1H','1H','1H'}; 
inter.zeeman.scalar={1.2 1.2 1.2 2.0};
inter.coupling.scalar=cell(4,4);
Jcell={0 -7.5 -7.5 0; -7.5 0 -7.5 0; -7.5 -7.5 0 0; 0 0 0 0};
inter.coupling.scalar=Jcell;
% set displacement of fouth spin
x = 3.0;
inter.coordinates = {[0.52 -0.90 0.36],[0.52 0.90 0.36],[-1.04 0.0 0.36],[x 0.0 0.36]};

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level = 4;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';

inter.rlx_keep='labframe';
inter.temperature=298;
inter.tau_c={6e-9};

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get thermal equilibrium state
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho_eq=equilibrium(spin_system,H);

% Build the relaxation superoperator
% R = relaxation(spin_system);
% or load a relaxation matrix and modify with thermalize
load('R_gce_10800k_6ns_5ps_lev_4_ave_3Tc','R_gce');
R = R_gce;
% Apply regularisation
if exist('reg','var')&&(reg~=0)
    R=R-reg*unit_oper(spin_system);
end
% Make sure the unit state is not damped
U=unit_state(spin_system);
R=R-(U'*R*U)*(U*U');

T = inter.temperature;
R = thermalize(spin_system,R,H,T,rho_eq,'dibari');

% Start in a state with one spin inverted (remote spin)
Lz4=state(spin_system,{'Lz'},{4});
rho=rho_eq-2*Lz4*(Lz4'*rho_eq)/norm(Lz4)^2;

% Compute the evolution trajectory
coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2}) ...
    state(spin_system,{'Lz'},{3}) state(spin_system,{'Lz'},{4})];
answer = evolution(spin_system,1i*R,coil,rho,0.05,100,'multichannel');


%% Plot the longitudinal magnetization
figure(); 
plot(linspace(0,5,101),answer/4.835e-7);
xlabel('time, seconds'); grid on;
ylabel('% equilirium intensity')
%end