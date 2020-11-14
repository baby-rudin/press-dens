clear;clc;

% parameters
m_O = 16;
m_H = 1;

fqq = 138.95;

q_O     = -0.834;
q_H     = 0.417;
epsla_O = 0.64;
epsla_H = 0.192;
sigma_O = 0.354;
sigma_H = 0.0449;

q_h2o     = [q_O q_H q_H];
epsla_h2o = [epsla_O epsla_H epsla_H];
sigma_h2o = [sigma_O sigma_H sigma_H];

dr  = 0.03;
ndr = 100;

nH2O    = 500;
nFrame  = 1000;


%VMD提取的1001-2000间结构
fileID = fopen('H2O500.pdb','r');
formatSpec = ' %*s %*s %*s %*s %*s %f %f %f %*s %*s %*s';
cont = textscan(fileID,formatSpec);
x = cont{1};
y = cont{2};
z = cont{3};

% remove NaN
x((0:nFrame)*(3*nH2O+1)+1) = [];
y((0:nFrame)*(3*nH2O+1)+1) = [];
z((0:nFrame)*(3*nH2O+1)+1) = [];

% remove redundant structure
x = x(1:3*nH2O*nFrame, :);
y = y(1:3*nH2O*nFrame, :);
z = z(1:3*nH2O*nFrame, :);

% change Unit of length: nm
x = x ./ 10.0;
y = y ./ 10.0;
z = z ./ 10.0;

x_1 = reshape(x, [3*nH2O nFrame]);
y_1 = reshape(y, [3*nH2O nFrame]);
z_1 = reshape(z, [3*nH2O nFrame]);

% atom weight matrix
mass_mask = zeros(3*nH2O, nFrame);
mass_mask((1:nH2O)*3 - 2,:) = m_O;
mass_mask((1:nH2O)*3 - 1,:) = m_H;
mass_mask((1:nH2O)*3 - 0,:) = m_H;

% calculate mass center
x_center = sum(mass_mask .* x_1) ./ sum(mass_mask);
y_center = sum(mass_mask .* y_1) ./ sum(mass_mask);
z_center = sum(mass_mask .* z_1) ./ sum(mass_mask);

% generate move vector
move_vecx = repmat(x_center, 3*nH2O, 1);
move_vecy = repmat(y_center, 3*nH2O, 1);
move_vecz = repmat(z_center, 3*nH2O, 1);

% move mass center of the system to origin
x_2 = x_1 - move_vecx;
y_2 = y_1 - move_vecy;
z_2 = z_1 - move_vecz;

% calculate mass center of each h2o
h2o_x_center = reshape(sum(reshape(x_2.*mass_mask, 3, nH2O, nFrame)), nH2O, nFrame) / (m_O + m_H *2);
h2o_y_center = reshape(sum(reshape(y_2.*mass_mask, 3, nH2O, nFrame)), nH2O, nFrame) / (m_O + m_H *2);
h2o_z_center = reshape(sum(reshape(z_2.*mass_mask, 3, nH2O, nFrame)), nH2O, nFrame) / (m_O + m_H *2);

% calculate distance between from mass center of each h2o to point O (mass center of system)
distance = sqrt(h2o_x_center.^2 + h2o_y_center.^2 + h2o_z_center.^2);

% calculate average density of molecule in each shell
nNum = zeros(1, ndr);
nDen = zeros(1, ndr);

for i= 1 : ndr
    inR = dr * (i-1);
    exR = dr * i;
    dV  = 4/3 * pi * (inR^3 - exR^3);
    match_idx = find(distance > inR & distance < exR);
    
    nNum(i) = length(match_idx);
    nDen(i) = nNum(i) / dV;
end

nNum = nNum / nFrame;
nDen = nDen / nFrame;


% calculate pressure U
pressU = zeros(ndr, nFrame);

parfor m = 1 : nFrame
    for i = 1 : nH2O - 1
    for j = i+1 : nH2O
        
        Pi = [x_2(3*i-2:3*i,m)  y_2(3*i-2:3*i,m)  z_2(3*i-2:3*i,m)];
        Pj = [x_2(3*j-2:3*j,m)  y_2(3*j-2:3*j,m)  z_2(3*j-2:3*j,m)];

        com_i = [h2o_x_center(i,m) h2o_y_center(i,m) h2o_z_center(i,m)];
        com_j = [h2o_x_center(j,m) h2o_y_center(j,m) h2o_z_center(j,m)];

        ri = distance(i,m);
        rj = distance(j,m);

        Pij = com_j - com_i;
        rij = sqrt(Pij*Pij');
        
        % Joined force
        fij = zeros(1,3);
            
        for ia = 1:3
        for ja = 1:3
            
            Pa = Pi(ia,:);
            Pb = Pj(ja,:);
                    
            sigmaij = (sigma_h20(ia) + sigma_h2o(ja)) / 2;
            epslaij = sqrt(epsla_h2o(ia) * epsla_h2o(ja));
                    
            Pab = Pb - Pa;
            rab = sqrt(Pab*Pab');
                    
            force = epslaij * 12 * (sigmaij^12 / rab^13 - sigmaij^6 / rab^7)...
                    + fqq * q_h2o(ia) * q_h2o(ja) / rab^2;
            
            fij = fij + force * Pab / rab;
                    
        end % end of ia
        end % end of ja
            
        % 总力在质心上的投影
        ffij = (fij * Pij') /rij;
            
        % distance from line (i--j) to mass center (origin)
        d = sqrt(ri^2 - (com_i * Pij')^2 / rij^2);
        
        % for each Spherical shell,
        % check whether ffij has influence on this shell
        for idr = 1 : ndr
            bR = idr * dr;
            cosine = (bR > d)*sqrt(bR^2 - d^2)/bR + (bR <= d)*0.0;
            
            if (ri - bR) * (rj - bR) < 0.0
                pressU(idr,m) = pressU(idr,m) + ffij * cosine * 1 / (4.0 * pi * bR^2);
            elseif ri > bR && rj > bR && abs(ri^2 - rj^2) < bR^2 && d < bR
                pressU(idr,m) = pressU(idr,m) + ffij * cosine * 2 / (4.0 * pi * bR^2);
            end
        end % end of idr
        
    end % end of loop i
    end % end of loop j
    
end % end of loop m (for each frame)

R   = (1:ndr)' * dr;
PU  = sum(pressU,2) / nFrame;
plot(R,PU)
