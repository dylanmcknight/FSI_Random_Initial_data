%This program inputs a square mesh and outputs the Taylor Hood mesh and
%Struct mesh
 
%Is working

% clear
% close all

%% Load in the square mesh
%  Square41
%    Square145
%    Square545
%     Square2113
   Square8321
%    Square33025

%% Define the relevant mesh parameters

Mf = msh.nbNod;
Lf = length(msh.TRIANGLES6(:,1));
Tf = length(msh.TRIANGLES6(1,:))-1;

J0f = msh.LINES3(:,1:3);
J0f = unique(J0f);
M0f = length(J0f);

Jf = transpose(setdiff(1:Mf,J0f));
Nf = length(Jf);

%% Set up the pressure mesh.

Jp = unique(msh.TRIANGLES6(:,1:3)); %Index set to translate between pressure nodes and fluid nodes
Mp = length(Jp);

msh.TH3 = zeros(Lf,4);

for ell = 1:1:Lf
    for r = 1:1:3
        msh.TH3(ell,r) = find(Jp == msh.TRIANGLES6(ell,r));
    end
end

%% Builds the Plate Mesh in terms of Fluid Mesh node numbers, then builds a local Plate Mesh
msh.STRUCT3 = [];

for i = 1:1:length(msh.LINES3(:,1))
    if msh.POS(msh.LINES3(i,3),2) == 1 %If the y component of the third node in the line segment is 1, then the line element must be on the top of the square
        msh.STRUCT3 = [msh.STRUCT3; msh.LINES3(i,:)]; %We will need to call the third node when computing the Dirichlet maps. Besides this it is unneeded
    end
end

Omega = unique(msh.STRUCT3(:,1:3)); %fluid node indicies on the plate portion of boundary
Omega_tilde = unique(msh.STRUCT3(:,1:2));
bdOmega = []; %Nodes on boundary of Omega where clamped boundary conditions are in effect



%This will define the portion of the boundary on which the problem has
%Dirichlet B.C.'s.
for i = 1:1:length(Omega_tilde)
    if (msh.POS(Omega_tilde(i),1) == 0) || (msh.POS(Omega_tilde(i),1) == 1) %The left and right endpoints
     bdOmega = [bdOmega, Omega_tilde(i)]; 
    end
end 

h = abs(msh.POS(msh.STRUCT3(1,1),1) - msh.POS(msh.STRUCT3(1,2),1));
Ls = length(msh.STRUCT3(:,1));

%% Makes the Local Plate Mesh
msh.LOCSTRUCT = zeros(Ls,5);

%This renumbers the nodes of msh.STRUCT3 to be 1:Mstruct
for i = 1:1:Ls
   for j = 1:1:2
     msh.LOCSTRUCT(i,j) = find(Omega_tilde == msh.STRUCT3(i,j));
   end
end

msh.LOCSTRUCT(:,3) = Ls+2:2:3*Ls+1;
msh.LOCSTRUCT(:,4) = msh.LOCSTRUCT(:,3) + 1;

Ms = 3*Ls+1; %Number of Plate Nodes
Ts = length(msh.LOCSTRUCT(1,:))-1; %Number of plate nodes per element
dofs = Ms + (Ms+2)/3; %Number of plate degrees of freedom.

M0s = length(bdOmega); 
J0s = []; %Renumbered Node indicies on the boundary of calO

for i = 1:1:M0s
   if msh.POS(Omega_tilde(i),1) == 0 || msh.POS(Omega_tilde(i),1) == 1
      J0s = [J0s, i]; 
   end
end

Js = setdiff(1:Ms,J0s); %Node indicies on the interior of Omega
Ns = length(Js);

%This finds the associated fluid triangular element for each linear
%structure element.

msh.LINE_TO_TRIANGLE = zeros(Ls,Tf+1);
for i = 1:1:Ls
    for j = 1:1:Lf
       if length(setdiff(msh.TRIANGLES6(j,:),msh.STRUCT3(i,:))) == 3
            msh.LINE_TO_TRIANGLE(i,:) = msh.TRIANGLES6(j,:);
       end
    end
end

%% Defines the positions of the Plate Nodes
msh.LOCPOS = zeros(Ms,2);

for ell = 1:1:Ls 
    x1 = msh.POS(msh.STRUCT3(ell,1),1);
    msh.LOCPOS(msh.LOCSTRUCT(ell,1),1) = x1;
    msh.LOCPOS(msh.LOCSTRUCT(ell,3),1) = x1 + h/3;
    msh.LOCPOS(msh.LOCSTRUCT(ell,4),1) = x1 + 2*h/3;
end

msh.LOCPOS(msh.LOCSTRUCT(Ls,2)) = msh.LOCPOS(msh.LOCSTRUCT(Ls,1)) + h;

clear x1

%% Defines the Hermite pts on in terms of both local and global node numbers.

msh.HERMITE = msh.LOCSTRUCT(:,1:2);
msh.HERMITE = unique(msh.HERMITE(:));


clear ell i j r