
%% Coordinate extract and R matrix calcualtion from an AMBER trajectory
%   J. Prestegard 01/09/20 03/17/20 04/08/20 07/15/20
%   Uses functions from mdtoolbox - Yasuhiro Matsunaga june 19, 1919
%   example is from a 10 ns trajectory on sucrose
%function Relaxation_matrix_ngce_from_trajectory_JHP_IK7()
%
%Give filename for EXCEL workbook used for data input and output and for Relaxation matrices.
% Workbook should have a sheet named 'Shifts' with chemical shifts in NMR
% Star format and a blank sheet named 'Output'
filename = 'Sucrose_opc_H1_6spin_2k-3kns_2ps_IK7_21T_24w';
% enter prmtop file name
prmtop_name = 'sucrose.opc.bare.prmtop';
% enter trajectory coordinate file name
crd_file = 'sucrose.opc.4md.2k-3k.ns.noWAT.crd';
% enter time step in seconds (time between saved frames)
t_step = 2e-12;
%enter frame increment (to skip some frames)
fr_inc = 1;
t_step = t_step*fr_inc;
% select residue and atom of interest
residue = 2;
atom = 'H1';
% select reference frame for calcualting distances
ref_frame = 255500;
ref_frame = ceil(ref_frame/fr_inc);
% specify cutoff distance
rcut = 4.0;
% adddional residues to be added
%residue_add = [1,2];
%atom_add = {'H3','H4'};
% specify number of workers
nwrk = 12;
%%
prmtop=readprmtop(prmtop_name);
%   get natom from prmtop output or use prmtop.natoms
%   get all coordinates excluding box 
[trj_all]=readmdcrd(prmtop.natom,crd_file);
% get number of frames
n_frames = size(trj_all,1);
if fr_inc > 1
    trj_all = trj_all(1:fr_inc:n_frames,:);
    n_frames = size(trj_all,1);
end
% write this to a .mat file
% need to append box matrix if box coordinates are at end of each frame.
save('trj_all_2','trj_all');
    
%% Get protons within a cutoff distance = of a selected atom
% Get coords of selected atom
index_resid = selectid(prmtop.residue_id,residue);
index_H = selectname(prmtop.atom_name, atom);
index_cent = index_resid & index_H;
cent_id = find(index_cent); % atom to be observed
% Get coordinates for observed atom from ref_frame
crd2 = [trj_all(ref_frame,cent_id*3-2), trj_all(ref_frame,cent_id*3-1), trj_all(ref_frame,cent_id*3)]; 
% Make file of H coordinates
index_all_H = selectname(prmtop.atom_name, 'H*');
% Eliminate H*O if sample in D2O
index_all_O = selectname(prmtop.atom_name, '*O');
index_all_OH = index_all_H & index_all_O;
index_all_H = index_all_H - index_all_OH;
all_H=find(index_all_H);
n = size(all_H);
% get coordinates for all non O*H H
H_coords = zeros(1,n(1,1)*3);
k=0;
for i = all_H.'
    for j = 1:3
        H_coords(k+j) = trj_all(ref_frame,(i-1)*3+j);
    end
    k=k+3;
end
%%
[pair, dist] = searchrange_exhaustive(H_coords, crd2, rcut); 

%%   get the coordinates for selected atom plus all within cutoff
pair_trn=pair.';
cut_list=pair_trn(2,:).';
np = size(cut_list);
ntp = size(trj_all,1);
%ntp = 5000; % limit number of frames if desired
if ntp > n_frames
    ntp = n_frames;
end
% turn cutlist into trajectory index calls
cut_list2 = zeros(1,np(1));
for i = 1:np(1)
    cut_list2(i) = all_H(cut_list(i,1));
end
%% add items to cut_list;
%size_add = size(residue_add,2);
%for i=1:size_add
    %index_res_add = selectid(prmtop.residue_id, residue_add(1,i));
    %index_atom_add = selectname(prmtop.atom_name, atom_add{i});
    %index_add = index_res_add & index_atom_add;
    %cut_list2(np(1)+i) = find(index_add);
%end
np = size(cut_list2,2);
% find coordinates in Spinach format
coords_cut = cell(ntp,np(1));
for ii = 1:ntp
    k=1;
    for i=cut_list2
        coords_cut{ii,k} = [trj_all(ii,(i-1)*3+1),trj_all(ii,(i-1)*3+2),trj_all(ii,(i-1)*3+3)]; 
        k=k+1;
    end
end
% save ref frame coordinates
coords_ref = coords_cut(ref_frame,:);
save('Coords_cut','coords_cut','-v7.3');
save('Cut_list','cut_list2');
% create a matfile variable to allow extraction of one trajectory segment at a time
trj_cut = matfile('Coords_cut');

%% Provide cut_list2 and atom to be observed to script that writes atom types, shifts and couplings 
% to EXCEL workbook to be called by Spinach simulation scripts
wb_filename = strcat(filename,'.xlsx');
sheetname_shifts = 'Shifts';
% Feed atom information from prmtop and cut_list; get shifts for checking
% Prepare_xlsx adds tries to fill matrix, but can pause to add couplings manually
zeeman = Prepare_xlsx_JHP_2b(wb_filename,sheetname_shifts,prmtop_name,cut_list2,trj_all(ref_frame,:),dist);
clear('trj_all');
% Read data from EXCEL work book for Spinach simulation
sheetname_output = 'Output';
data = readcell(wb_filename,'Sheet',sheetname_output);

%% Spinach part begins

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0'; %'none' 
bas.level=3;

% System specification
sys.magnet=21.2;
sys.isotopes=cell(1,np(1));
sys.isotopes(1:end)={'1H'};

% Chemical shifts
inter.zeeman.scalar=data(2:(np(1)+1),5)';

% Scalar couplings
inter.coupling.scalar=data(2:(np(1)+1),13:(13+np-1));
%inter.coupling.matrix{np(1),np(1)}=zeros(3);

% Non-MD relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='zero';
inter.rlx_keep='labframe';
inter.tau_c={1.0e-10};
inter.temperature=298;

% Run Spinach housekeeping
spin_system=create(sys,inter);

spin_system=basis(spin_system,bas);

% Set laboratory frame
spin_system=assume(spin_system,'labframe');

% Calculate H0
H0=hamiltonian(spin_system);

% Get size of trajectory and define segment size
trj_size = size(coords_cut,1);

% Trajectory length in s
trj_length = t_step*trj_size;

% Estimated rotational correlation time
TauC = 1.0e-10;

% Find # frames in 3*TauC
seg_size = floor(trj_size*3*TauC/trj_length); 
nsegs = floor(trj_size/seg_size);

% Preallocate the answer
R_gce=zeros(size(H0));

% Hush up Spinach
spin_system.sys.output='hush';
disp('Spinach output hushed');

% Loop over segments
% set number of workers
delete(gcp('nocreate'));
parpool(nwrk);
parfor nn=1:nsegs
    
    % Grab a local copy of the data structure
    ss_local=spin_system;
    
    % Preallocate segment array
    H1_seg = cell(seg_size,1);
    
    % Get a trajectory segment
    trj_seg = trj_cut.coords_cut((nn-1)*seg_size+1:nn*seg_size,:); %#ok<PFBNS>
    
        % Loop over the segment
        for n=1:seg_size
            
            % Set current coordinates
            ss_local.inter.coordinates=trj_seg(n,:);
            ss_local.inter.pbc={};
            
            % Ingest the coordinates
            ss_local_frame=dipolar(ss_local);
            
            %  Get H1 for each frame in segment
            [~,Q]=hamiltonian(ss_local_frame);
            H1_seg{n}=orientation(Q,[0.0,0.0,0.0]);
        
        end
        
    % Get the GCE relaxation matrix
    R_gce=R_gce+ngce_IK_JHP_3(ss_local,H0,H1_seg,seg_size,t_step);
    disp(['finished computing segment: ' num2str(nn)]); drawnow();
    
end

% Average and tidy up
R_gce=R_gce/nsegs;
R_gce=real(R_gce); 
R_gce=(R_gce+R_gce')/2;

% Get coords from reference frame of trajectory
spin_system.inter.coordinates=trj_cut.coords_cut(ref_frame,:);

spin_system.inter.pbc={};

% Ingest coordinates
spin_system=dipolar(spin_system);

% Get Redfield matrix for comparison
R_red=relaxation(spin_system);

% Save the files
save(strcat(filename,'_R_gce'),'R_gce');
save(strcat(filename,'_R_red'),'R_red');

% Get the answers to the user
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
disp(['RMSD Numerical vs Analytical: ' num2str(mean(abs(nonzeros(R_gce-R_red)))) ' Hz']);

%end

