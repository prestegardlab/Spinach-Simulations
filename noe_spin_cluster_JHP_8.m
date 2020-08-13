% Nuclear Overhauser effect on the central spin pair of a cluster defined
% by a cutoff distance. The spript depends on loading data files produced
% by additional scripts: A workbook with complete spin system description
% and an relaxation matrix simulated externally. This version inverts the
% spin at the center of the cluster.
% This is based on an initial 2-spin NOE script by: 
% i.kuprov@soton.ac.uk
% Modified by J. Prestegard 12/31/19 03/17/20 06/25/20 08/07/20

%%function noe_spin_cluster_JHP_7()
% Identify EXCEL work book and sheet containing info for Spinach simulation
filename = 'Sucrose_opc_H1_8spin_2k-3kns_2ps_IK7_12w_auto.xlsx';
sheetname = 'Output';
% Name of Relaxation matrix to be used - and name of variable
R_name = 'Sucrose_opc_H1_8spin_2k-3kns_2ps_IK7_12w_R_gce.mat';
Var_name = 'R_gce';
% Enter mixing times in seconds as a vector
mxt= [0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8]; 

%% read in data from workbook
data = readcell(filename,'Sheet',sheetname);
% Find size of spin system
% find number of comment lines
n = 0;
while ~isnumeric(data{n+1,1})
    n = n+1;
end
% find number of spins
n_spins = size(data,1)-n;

% Assign the spin type
isotope ={'1H'};
for j = 1:n_spins
    sys.isotopes(1,j)=isotope; 
end
inter.zeeman.scalar=data(n+1:(n+n_spins),5)';
inter.coupling.scalar=data(n+1:(n+n_spins),13:(12+n_spins));
for i=1:n_spins
    inter.coordinates{i}=cell2mat(data(n+i,9:11));
end

%% Set the spin system

% Magnet field
sys.magnet=14.1;

% Basis set
bas.formalism='sphten-liouv';
bas.approximation='IK-0';
bas.level=3;

% Relaxation theory parameters
inter.relaxation={'redfield'};
inter.equilibrium='dibari';
inter.rlx_keep='kite';
inter.temperature=298;
inter.tau_c={0.1e-9};

% Proximity cut-off
sys.tols.prox_cutoff=4.0;

% Spinach housekeeping
spin_system=create_JHP(sys,inter);
spin_system=basis(spin_system,bas);
H = hamiltonian(assume(spin_system,'labframe'),'left');
% Get thermal equilibrium state
rho_eq=equilibrium(spin_system,H);
%% Calculate NOE based on relaxation matrix
load(R_name,Var_name);
R=thermalize_JHP(spin_system,R_gce,H,298,rho_eq,'dibari');
% Invert non-selected spins one at a time and compute Z magnetization at dwells time dw
answer_f=zeros(1,n_spins);
% find inverted spin
for j = 2:n_spins+1
    if cell2mat(data(j,7))==0.0
        inv_spin=j-1;
    end
end
%% Prepare workbook for output
alpha=' A B C D E F G H I J K L M N O P Q R S T U V W X Y ZAAABACADAEAFAGAHAIAJAKALAMANAOAPAQARASATAUAVAWAXAYAZ';
data2 = readcell(filename,'Sheet',sheetname);
end_of_data=2*size(data2,2)+1;
col=extractBetween(alpha,end_of_data,end_of_data+1);
col_id=char(strcat(col,'1'));
NOE_Label = 'NOEmix s';
writematrix(NOE_Label,filename,'Sheet','Output','Range',col_id);
col=extractBetween(alpha,end_of_data+2,end_of_data+3);
col_id=char(strcat(col,'1'));
writematrix(mxt(:)',filename,'Sheet','Output','Range',col_id);
k = 0;
%%
for kmix = mxt
    k = k+1;
    for j = 2:n_spins+1
        ob_spin=j-1;
        Lz1=state(spin_system,{'Lz'},{inv_spin});
        rho=rho_eq-2*Lz1*(Lz1'*rho_eq)/norm(Lz1)^2;
        coil=[state(spin_system,{'Lz'},{ob_spin})]; % state(spin_system,{'Lz'},{inv_spin})];
        answer=evolution(spin_system,1i*R,coil,rho,kmix,1,'observable');
        if ob_spin == inv_spin
            answer_f(1,j-1)= -answer(2,1)/answer(1,1);
        else
            answer_f(1,j-1)=(answer(2,1)-answer(1,1))/answer(1,1);
        end
    end
%% write NOEs to EXCEL sheet
col=extractBetween(alpha,end_of_data+2*k,end_of_data+2*k+1);
col_id=char(strcat(col,'2'));
writematrix(answer_f',filename,'Sheet','Output','Range',col_id);
end
%end

