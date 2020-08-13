% Nuclear overhauser effect for a methyl plus one proton modeled with
% static pairs of pseudo spins. Assumptions are that the like spin states
% of a methyl group (a,a,a) can be modeled with spins having a normal gamma,
% and unlike spins (a,a,b) can be modeled with spins having a reduced gamma
% (~1/3). These states have relative abundances of 1:3 as refelcted in the
% outer and inner lines of a 1H coupled 13C detected multiplet of a methy.
%
% If one wishes to model a system with rapid averaging among the positions
% of methyl proton, one would simply do a 1/r^3 average over the three sites
% for each proton.  This works in this geometry since theta doesn't change. 
% This is also true for the on axis models and the in plane model with the 
% remote proton rotated by 90 degrees.
% For this model the scaling factor (comparing (ave(1/r^3))^2 to Sum(1/r^6)
% is 0.86.  for the on axis model it is 1.0, for the 90 rotated model 0.64
% Note that all scaling factors go to 1 at large distances to remote spin.
%
% Based on a 2_spin_NOE script by I. Kuprove i.kuprov@soton.ac.uk
% modified from noe_two_spin by J Prestegard 02/27/20
%
% 08/03/20 - changed to sum pairwise interations over all 6 contributions 

%function noe_static_spin_pair_Methyl_model()

% Set the spin system
sys.isotopes={'1H','1H','1H','1H'}; 
inter.zeeman.scalar={1.2 1.2 1.2 2.0};
inter.coupling.scalar=cell(4,4);
Jcell={0 -7.5 -7.5 0; -7.5 0 -7.5 0; -7.5 -7.5 0 0; 0 0 0 0};
inter.coupling.scalar=Jcell;
x=3.0; % variable postion of non-methyl proton
y=10.0; %displacement to minimize interactions
% need three spins positioned to sample the two possible methyl-remote spin
% separations with a second methyl spin at a distance and a methyl-methyl 
% separation with the remote spin at a distance.  This will allow only spin pair
% contributions to relaxation. The following should work: one methyl-methly pair at a separation
% of 1.8A and the remote spin at a distance, one methyl-remote pair at a separation
% of ((x-0.52)^2+0.91^2)^1/2 and a methyl spin at a distance, and one with a 
% methyl-remote pair at a separation of (x+1.04)^2 and a methyl spin at a distance 
coords = {[0.52+y -0.90 0.36],[0.52 0.90+y 0.36],[-1.04 0.0 0.36],[x 0.0 0.36];...
    [0.52 -0.90+y 0.36],[0.52 0.90 0.36],[-1.04+y 0.0 0.36],[x 0.0 0.36];...
    [0.52 -0.90+y 0.36],[0.52 0.90 0.36],[-1.04 0.0 0.36],[x+y 0.0 0.36];...
    [0.52+y -0.90 0.36],[0.52 0.90+y 0.36],[-1.04+y 0.0 0.36],[x 0.0 0.36];...
    [0.52+y -0.90 0.36],[0.52 0.90+y 0.36],[-1.04 0.0 0.36],[x+y 0.0 0.36];...
    [0.52+y -0.90 0.36],[0.52 0.90 0.36],[-1.04+y 0.0 0.36],[x+y 0.0 0.36]};
nspins = size(coords,2);
nstates= size(coords,1);
R_size=2^(2*nspins);
% weights to be used in summing relaxation contributions
% scale factor for gamma in rotationally averaged model
gamma_scale=0.333;
% Initialize R matricies
R_like_sum = zeros(R_size);
R_unlike_sum = zeros(R_size);

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='none';
% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';

inter.rlx_keep='labframe';
inter.temperature=298;
inter.tau_c={6e-9};

for i=1:nstates
inter.coordinates=coords(i,1:nspins);

% Spinach housekeeping
spin_system=create(sys,inter);
spin_system=basis(spin_system,bas);

% Get thermal equilibrium state
H=hamiltonian(assume(spin_system,'labframe'),'left');
rho_eq=equilibrium(spin_system,H);

% Build the relaxation superoperator
R_temp = relaxation(spin_system);
R_like_sum = R_like_sum + R_temp;

% Scale relaxation by gamma effective ^2 for inner line
R_unlike_sum = R_unlike_sum + R_temp*(gamma_scale)^2;
end
% Optionally weight by factor to account for 1/r^3 averaging if non-static
%R_like_sum = R_like_sum .* 0.86;
%R_unlike_sum = R_unlike_sum .* 0.86;
% Start in a state with one spin inverted (remote spin)
Lz2=state(spin_system,{'Lz'},{4});
rho=rho_eq-2*Lz2*(Lz2'*rho_eq)/norm(Lz2)^2;

% Compute the evolution trajectory
coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2}) state(spin_system,{'Lz'},{3})...
    state(spin_system,{'Lz'},{4})];
answer_in=evolution(spin_system,1i*R_unlike_sum,coil,rho,0.05,100,'multichannel');
answer_out=evolution(spin_system,1i*R_like_sum,coil,rho,0.05,100,'multichannel');
% weight buildup curves by degeneracies unlike and outer like spin states.
answer=0.75*answer_in+0.25*answer_out;
%% Plot the longitudinal magnetization
% plot shows each methyl proton individulally - must add to get total
% methyl intensity.
figure(); 
plot(linspace(0,5,101),answer/(4.835e-7));
xlabel('time, seconds'); grid on;
ylabel('% equilirium intensity')
hold off
%end

