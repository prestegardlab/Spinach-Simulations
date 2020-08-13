% This script reads chemical NMR Star format shifts from an existing EXCEL sheet
% First lines can be comments; remainder will be lines that begin with atom numbers
% followed by residue #, residue name,atom name (twice) and chemical shift.
% The script then prepares an "Output" sheet to store shifts and addtional data used
% by or returned from Spinach simulatiions. A scalar coupling constant
% matrix is written using extracted dihedral angles and a Karplus relation
% to get 3-bond couplings and inserting a default -15 Hz for 2-bond couplings.
% Test case is sucrose
%
% J. Prestegard, 01/14/20, 1/29/20, 04/23/20
%

function [zeeman] = Prepare_xlsx_JHP_2b(filename,sheetname,prmtop_name,cut_list,coords,dist)
%% Read shift data from EXCEL work book for Spinach simulation
% filename ans sheetname refer to the EXCEL workbook used for data storage,
% prmtop_name is for the trajectory used, cutlist contains atom numbers for
% the cluster of protons under study, coords contians their coordinates and
% dist is the distance between the selected protons (center of cluster) and
% other spins
shifts = readcell(filename,'Sheet',sheetname);
% find number of comment lines
n = 1;
while ~isnumeric(shifts{n,1})
    n = n+1;
end
% find number of atoms in shift list
natm = 1;
while isnumeric(shifts{(n+natm),1})
    natm = natm+1;
end
% Get atom information from prmtop of trajectory and cut_list produced by
% Coordinate_extract script
prmtop=readprmtop(prmtop_name);
cutsize = size(cut_list,2);
zeeman = cell(1,cutsize);
res_num = zeros(1,cutsize);
atom_name = cell(1,cutsize);
res_name = cell(1,cutsize);
atom_id = zeros(1,cutsize);
% make residue number cells into 1 x natm matrix
res_num_all(n:natm) = cell2mat(shifts(n:natm,2));
nj = 0;
for j = cut_list  
    % extract residue id , residue name and atom_name
    jj = 0;
    nj = nj+1;
    res_num(nj) = prmtop.residue_id(j);
    atom_id(nj) = prmtop.atom_id(j);
    res_name{1,nj} = prmtop.residue_name(j,1:3);
    atom = prmtop.atom_name(j,1:4);
    for k = n:natm
        if res_num(nj) == res_num_all(1,k)
            atom_all = pad(cell2mat(shifts(k,4)),4);
            % try matching just 3 for PPM1 output
            if strncmpi(atom,atom_all,3)
                zeeman{1,nj} = shifts{k,6};
                jj = 3;
            end
        end
    end
        % if three character match is not found, try 2 characters
    if jj == 0
        for kk = n:natm       
            if res_num(nj) == res_num_all(1,kk)
                 atom_all = pad(cell2mat(shifts(k,4)),4);
                 if strncmpi(atom,atom_all,2)
                     zeeman{1,nj} = shifts{k,6};
                     jj = 2;
                 end
            end
        end
    end 
    atom_name{nj} = atom(1:4);
end
res_name = char(res_name);
%save('cut_list_zeeman','zeeman');
%save('cut_list_res_num','res_num');
%save('cut_list_atom_name','atom_name');
%save('cut_list_res_name','res_name');
writematrix(atom_id.',filename,'Sheet','Output','Range','A2');
writematrix(res_num.',filename,'Sheet','Output','Range','B2');
writematrix(res_name,filename,'Sheet','Output','Range','C2');
atom_typs = cell(1,cutsize);
for i = 1:cutsize
    atom_typs{i} = '1H';
end
writecell(atom_typs.',filename,'Sheet','Output','Range','D2');
writecell(zeeman.',filename,'Sheet','Output','Range','E2');
writecell(atom_name.',filename,'Sheet','Output','Range','F2');
writematrix('coupling',filename,'Sheet','Output','Range','M1');
writematrix(dist',filename,'Sheet','Output','Range','G2');
writematrix('chem shift',filename,'Sheet','Output','Range','E1');
writematrix('distance',filename,'Sheet','Output','Range','G1');
%get coordinates in reference frame for cutlist and write to workbook
writematrix('coordinates',filename,'Sheet','Output','Range','I1');
k = 0;
for i = cut_list
    k = k+1;
    write_pt = strcat('I',int2str(1+k));
    writematrix(coords(i*3-2:i*3), filename, 'Sheet', 'Output', 'Range', write_pt);
end
%% Make coupling matrix for Spinach simulation
J = zeros(cutsize);
% Get list of dihedrals from prmtop
dihedrals = prmtop.dihedrals_inc_hydrogen;
for j = 1:cutsize
% Convert cutlist id to dihedral pointer
dihed_id = (cut_list(j)-1)*3;
dihed_size = size(dihedrals,1);
% Search for dihed_id in first dihedral position and a cut_list proton in fouth
for i=1:dihed_size
    if dihedrals(i,1)==dihed_id
        for ii = 1:cutsize
            if dihedrals(i,4)==(cut_list(ii)-1)*3
                tort = dihedrals(i,1:4);
                tort = tort/3+1;
                dihed_angle = calcdihedral(coords,tort(1:4));
                J(j,ii) = 8.14*cos(dihed_angle)^2 - 0.61*cos(dihed_angle) - 0.15; 
            end
        end
    end
end
end
%% Get list of angles from prmtop
angles = prmtop.angles_inc_hydrogen;
for j = 1:cutsize
% Convert cutlist id to dihedral pointer
angle_id = (cut_list(j)-1)*3;
angle_size = size(angles,1);
% Search for angle_id in first angle position and a cut_list proton in third
for i=1:angle_size
    if angles(i,1)==angle_id
        for ii = 1:cutsize
            if angles(i,3)==(cut_list(ii)-1)*3
                J(j,ii) = -15;
            end
        end
    end
end
end  
%% Write coupling matrix to EXCELL spreadsheet'
writematrix(J,filename,'Sheet','Output','Range','M2');

end
