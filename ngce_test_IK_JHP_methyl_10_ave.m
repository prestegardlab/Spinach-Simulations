% Numerical integral route to a relaxation matrix (R_gce).
% The program included the option to average Hamiltonians to mimic fast
% motion between integral steps.
% Results are plotted against an analytical Redfileld matrix (R_red).
% Both matricies are stored
%
% Code is based on an intial version of ngce_test developed by I. Kuprov
% (i.kuprov@soton.ac.uk)
% modified by J. Prestegard 11/22/19 and tested 02/09/20 to average rotamer states of a methyl
% modified by J. Prestegard 02/24/20 to call ngce subroutine ncge_IK_JHP_3 by segment.

%%function ngce_test_IK_JHP_methyl_10_ave()
rng('shuffle');
% specify number of workers
nwrk = 24;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='secular';
inter.tau_c={6.0e-9};
inter.temperature=298;

%% Get Euler angles for a random walk on a sphere
tc_steps=1200; % Number of time steps in one tau_c
tau_c=inter.tau_c{1}; dt=tau_c/tc_steps;
traj_pts=tc_steps*3; % expected duration of correlation
npoints = 10800000; % number of points in trajectory must be a multiple of traj_pts
eulers=rwalk_IK(npoints,tau_c,dt); 
% break eulers into segments
% nsegs must be greater than 200
nsegs = npoints/traj_pts;
eulers_seg=zeros(traj_pts,3);

%% set possible coordinates 
coords = {[0.52 -0.90 0.36], [0.52 0.90 0.36], [-1.04 0.0 0.36], [3.0 0.0 0.36]; 
[0.52 0.90 0.36], [-1.04 0.0 0.36], [0.52 -0.90 0.36],[3.0 0.0 0.36];
[-1.04 0.0 0.36], [0.52 -0.90 0.36], [0.52 0.90 0.36], [3.0 0.0 0.36]};
num_states=size(coords,1);
% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0'; %'none' 
bas.level=4;
% System specification
sys.magnet=14.1;
sys.isotopes={'1H','1H','1H','1H'};
nspins=size(sys.isotopes,2);
inter.zeeman.scalar={1.2 1.2 1.2 2.0};
inter.coupling.scalar={0 -7.5 -7.5 0 ; -7.5 0 -7.5 0 ; -7.5 -7.5 0 0 ; 0 0 0 0};
% get initial Hamiltonian
inter.coordinates=coords(1,1:nspins,1);
spin_system=create_JHP(sys,inter);
spin_system=basis(spin_system,bas);
[H0,Q]=hamiltonian_JHP(assume(spin_system,'labframe'));
% reserve Hamiltonian and R space
H1i=cell(num_states,traj_pts,1);
H1_seg=cell(1,traj_pts);
R_gce=zeros(size(H0));
% set number of workers
delete(gcp('nocreate'));
parpool(nwrk);
% Hush up Spinach
spin_system.sys.output='hush';
disp('Spinach output hushed');
%% Use Euler angles segment by segment
for nn=1:nsegs
    first_pt=(nn-1)*traj_pts+1;
    last_pt=nn*traj_pts;
    eulers_seg(1:traj_pts,1:3)=eulers(first_pt:last_pt,1:3);
    %  Get Hamiltonian and Q for each rotmer
    for i=1:num_states
        inter.coordinates=coords(i,1:nspins,1);
        spin_system=create_JHP(sys,inter);
        spin_system=basis(spin_system,bas);
        [H0,Q]=hamiltonian_JHP(assume(spin_system,'labframe'));
        % rotate Q to each segment of Euler points
        parfor n=1:traj_pts   
            H1i{i,n} = orientation(Q,eulers_seg(n,:));
        end
    end
    % average H1
    for j=1:traj_pts
        H1_sum = zeros(size(H0));
        for k=1:num_states
            H1_temp = H1i(k,j);
            H1_sum = H1_sum + H1_temp(1,:);
        end
        H1_seg{1,j} = cell2mat(H1_sum)/num_states;
    end
    disp(['finished computing segment: ' num2str(nn)]); drawnow();
    %report(spin_system,' ');
    %report(spin_system,['finished computing segment: ' num2str(nn)]);
    
%% Get the GCE relaxation matrix
R_gce= R_gce+ngce_IK_JHP_3(spin_system,H0,H1_seg,traj_pts,dt);

end
% Average and Tidy up relaxation matrix
R_gce=R_gce/nsegs;
R_gce=real(R_gce); R_gce=(R_gce+R_gce')/2;
R_red=relaxation(spin_system);
save('R_red_10800k_6ns_5ps_lev_4_ave_3Tc_secular','R_red');
save('R_gce_10800k_6ns_5ps_lev_4_ave_3Tc_secular','R_gce');
%% Get the answers to the user
disp('Relaxation superoperator, numerical'); disp(full(R_gce));
disp('Relaxation superoperator, analytical'); disp(full(R_red));
% Do diagnostic plotting
figure();
gce_rates=diag(R_gce);
red_rates=diag(R_red);
min_rate=min([R_gce; R_red]);
max_rate=max([R_gce; R_red]);
plot(gce_rates,red_rates,'ro'); hold on;
plot([min_rate max_rate],[min_rate max_rate],'b-');
axis tight; box on; grid on;
xlabel('Numerical relaxation rates');
ylabel('Analytical relaxation rates');

% Get an accuracy indicator
disp(['RMSD vs Spinach: ' num2str(mean(abs(nonzeros(R_gce-R_red)))) ' Hz']);

%end

% stop the parallel pool
delete(gcp('nocreate'));
%end
