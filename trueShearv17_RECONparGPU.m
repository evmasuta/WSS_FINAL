function trueShearv17_RECONparGPU(sourceDir,outDir,expt,dimT,numFiles,parallelitySwitch,parallelityCutoff,onlyVel)

%% PART1 GEOMETRY INITIALIZATION
% Initialization variables
topDir = pwd;
addpath(genpath(topDir));
numFiles = double(str2num(numFiles));
startT = 1;
dimT = double(str2num(dimT));
startZ = 1;
dimZ = double(numFiles/dimT);

parallelitySwitch = double(str2num(parallelitySwitch));
parallelityCutoff = double(str2num(parallelityCutoff));
onlyVel = double(str2num(onlyVel));

% viscosity value taken from 10.1002/jmri.24959 (van OOij et al., 2015)
% Units of Pa*s == kg / (m * sec)
% At its top level, let's try to force everything into centimeter units
viscosity = 3.2*10^-3;
% convert to cPa
viscosity = viscosity * 100;

% COORDINATE PRE-PROCESSING: Use Intensity files
name_prefix = strcat(sourceDir,'/mag/Cph_');
name_suffix = '.mag';

firstFrameName = sprintf('%s%03d%s%03d%s',name_prefix,startT-1,'_Sec_',startZ,name_suffix);
lastFrameName = sprintf('%s%03d%s%03d%s',name_prefix,startT-1,'_Sec_',dimZ,name_suffix);

% Stores coordinate data for first and last slices at a given time
firstFrame = dicominfo(firstFrameName);
lastFrame = dicominfo(lastFrameName);
rowCol = dicomread(firstFrame);

% Angular conversion row/col/slice to x,y,z
temp = double(firstFrame.ImageOrientationPatient);
F(:,1) = temp(4:6);
F(:,2) = temp(1:3);
clear temp

% (x,y,z) of first voxels
T1 = double(firstFrame.ImagePositionPatient);
Tn = double(lastFrame.ImagePositionPatient);
n = dimZ;

% 3D affine constants
k1 = (T1(1) - Tn(1)) / (1-n);
k2 = (T1(2) - Tn(2)) / (1-n);
k3 = (T1(3) - Tn(3)) / (1-n);

% Row/col pixel spacing: dr is vertical spacing; dc is horizontal spacing
dr = double(firstFrame.PixelSpacing(1));
dc = double(firstFrame.PixelSpacing(2));
ds = double(firstFrame.SliceThickness);


% Affine matrix
A = double([F(1,1)*dr, F(1,2)*dc, k1, T1(1); F(2,1)*dr, F(2,2)*dc, k2, T1(2); F(3,1)*dr, F(3,2)*dc, k3, T1(3); 0, 0, 0, 1]);
%% Geometry storage: geom(:,:,:,1) = x, geom(:,:,:,2) = y, geom(:,:,:,3) = z
geom = zeros([size(rowCol),dimZ,3]);
temp = size(rowCol);
rowCol = temp;
for row = 1:1:temp(1)
    for col = 1:1:temp(2)
        for slice = 1:1:dimZ
            % Zero-indexed!
            rcs = double([row-1;col-1;slice-1;1]);
            coords = A*rcs;
            geom(row,col,slice,1) = coords(1);
            geom(row,col,slice,2) = coords(2);
            geom(row,col,slice,3) = coords(3);
        end
    end
end

[rubbish,dr] = max(abs(A(:,1)));
[rubbish,dc] = max(abs(A(:,2)));
[rubbish,ds] = max(abs(A(:,3)));
% NEED TO CONVERT TO cm TO MATCH UNITS OF VELOCITY (cm/s)
dr = A(dr,1) / 10;
dc = A(dc,2) / 10;
ds = A(ds,3) / 10;
spacings = [dr,dc,ds];

clear temp;

%% PART 2 INITIALIZE KERNELS
% tolerance value for rcs determination
tolerance = 0.001;
% Determination of rcs to x,y,z (intuitive)
dimDiff(1,:) = squeeze(geom(2,1,1,:) - geom(1,1,1,:));
dimDiff(2,:) = squeeze(geom(1,2,1,:) - geom(1,1,1,:));
dimDiff(3,:) = squeeze(geom(1,1,2,:) - geom(1,1,1,:));
for i = 1:1:3
    for j = 1:1:3
        if abs(dimDiff(i,j)) < tolerance
            dimDiff(i,j) = 0;
        end
    end
end

kernelsNorm = KernelMakerv5_Central_RECON(dimDiff, spacings);

xKernel = kernelsNorm{1};
yKernel = kernelsNorm{2};
zKernel = kernelsNorm{3};
delta = kernelsNorm{4};
clear kernelsNorm
% CHANGE KERNELS
kernelsVel = KernelMakerv4_RECON(dimDiff, spacings);
xKernelFwd = kernelsVel{1};
yKernelFwd = kernelsVel{2};
zKernelFwd = kernelsVel{3};
xKernelRvr = kernelsVel{4};
yKernelRvr = kernelsVel{5};
zKernelRvr = kernelsVel{6};
% should be redundant
%delta = kernelsVel{7};
clear kernelsVel
%% PART 3 file I/O, aka part that needs to be modded heavily for the cluster DONE

% geom prefix suffix
geom_prefix = strcat(sourceDir,'/mag/Cph_');
geom_suffix = '.mag';
% v0 prefix suffix
v0_prefix = strcat(sourceDir,'/rl/Cph_');
v0_suffix = '.rl';
% v1 prefix suffix
v1_prefix = strcat(sourceDir,'/ap/Cph_');
v1_suffix = '.ap';
% v2 prefix suffix
v2_prefix = strcat(sourceDir,'/si/Cph_');
v2_suffix = '.si';
% We loop by time to reduce the memory (RAM) load, which will hopefully
% prevent crashing and speed things up a bit...
parfor t = startT:1:dimT
    %clearing stuff here is key, otehrwise type of matrix could be carried
    %over (switch from double to int16)
    % initialize storage; no need for time dim
    allEdge = double(zeros([rowCol(1), rowCol(2), dimZ]));
    %     xShear = double(zeros([rowCol(1), rowCol(2), dimZ]));
    %     yShear = double(zeros([rowCol(1), rowCol(2), dimZ]));
    %     zShear = double(zeros([rowCol(1), rowCol(2), dimZ]));
    % Initialize intensity storage
    intensity = double(zeros([rowCol(1), rowCol(2), dimZ]));
    % Initialize velocity storage
    velocity = double(zeros([rowCol(1), rowCol(2), dimZ, 3]));
    % Initialize intensity gradient storage
    gradIntensity = double(zeros([rowCol(1), rowCol(2), dimZ]));
    % initialize velocity gradient storage
    gradVelocityEdge = double(zeros([rowCol(1), rowCol(2), dimZ, 3, 3, 2]));
    gradVelocity = double(zeros([rowCol(1), rowCol(2), dimZ, 3, 3]));
    % initialize information storage: data are: [intensity,v0,v1,v2]
    information = cell(dimZ,4);
    
    % Start reading in data
    for i=startZ:1:dimZ
        % Filenames
        geom_name = sprintf('%s%03d%s%03d%s',geom_prefix,t-1,'_Sec_',i,geom_suffix);
        v0_name = sprintf('%s%03d%s%03d%s',v0_prefix,t-1,'_Sec_',i,v0_suffix);
        v1_name = sprintf('%s%03d%s%03d%s',v1_prefix,t-1,'_Sec_',i,v1_suffix);
        v2_name = sprintf('%s%03d%s%03d%s',v2_prefix,t-1,'_Sec_',i,v2_suffix);
        % Read'n'store; N.B. that information for intensity and velocity
        % are stored in the same call array!!!
        information{i,1} = dicominfo(geom_name);
        intensity(:,:,i) = dicomread(geom_name);
        information{i,2} = dicominfo(v0_name);
        velocity(:,:,i,1) = dicomread(v0_name);
        information{i,3} = dicominfo(v1_name);
        velocity(:,:,i,2) = dicomread(v1_name);
        information{i,4} = dicominfo(v2_name);
        velocity(:,:,i,3) = dicomread(v2_name);
    end
    intensity = double(intensity); velocity = double(velocity);
    %% PART 4: CALCULATE GRAD INTENSITY (positive grad, inward normal since lumen is bright CHECK THIS)
    %d/dx
    speed = squeeze((velocity(:,:,:,1).^2 + velocity(:,:,:,2).^2 + velocity(:,:,:,3).^2).^0.5);
    gradIntensity(:,:,:,1) = convn(speed(:,:,:),xKernel,'same') / delta(1);
    % d/dy
    gradIntensity(:,:,:,2) = convn(speed(:,:,:),yKernel,'same') / delta(2);
    % d/dz
    gradIntensity(:,:,:,3) = convn(speed(:,:,:),zKernel,'same') / delta(3);
    % magnitude of grad Intensity
    magGradIntensity = (gradIntensity(:,:,:,1).^2 + gradIntensity(:,:,:,2).^2 + gradIntensity(:,:,:,3).^2).^0.5;
    magGradIntensity = repmat(magGradIntensity,1,1,1,3);
    magGradIntensity = double(magGradIntensity);
    % normal vector
    normal = gradIntensity ./ magGradIntensity;
    alpha = normal(:,:,:,1);
    beta = normal(:,:,:,2);
    gamma = normal(:,:,:,3);
    % Hadamard product and scaling; may change the way we scale later...
    scale_factor = 5000 * 100 * 2 / 4;
    allEdge(:,:,:) = double(intensity .* magGradIntensity(:,:,:,1));
    maxIntensity = double(max(max(max(max(allEdge)))));
    allEdge = double(allEdge / maxIntensity * scale_factor);
    
    %% PART 5 CALCULATE PARALLEL AND ORTHOGONAL COMPONENTS OF VELOCITY
    % Note the nomenclature here is the direction of the velocity relative
    % to the surface; normVel is pointing in the direction of the normal
    % while parVel is orthogonal to normal.
    % Orthogonal component; normal is a unit vector so no need to divide
    normVelMag = repmat(dot(velocity,normal,4),1,1,1,3);
    normVel = normVelMag .* normal;
    parVel = velocity - normVel;
    velocity = parVel;
    %% PART 5.1: NEW
    % Calculate speed
    speed_squared = velocity(:,:,:,1).^2 + velocity(:,:,:,2).^2 + velocity(:,:,:,3).^2;
    gradSpeed = double(zeros([rowCol(1), rowCol(2), dimZ, 3]));
    gradSpeed(:,:,:,1) = convn(speed_squared(:,:,:,1),xKernel,'same') / delta(1);
    gradSpeed(:,:,:,2) = convn(speed_squared(:,:,:,1),yKernel,'same') / delta(2);
    gradSpeed(:,:,:,3) = convn(speed_squared(:,:,:,1),zKernel,'same') / delta(3);
    % pseudo-angle from dot product
    normNorm = (normal(:,:,:,1).^2 + normal(:,:,:,2).^2 + normal(:,:,:,3).^2).^0.5;
    normSpeed = (gradSpeed(:,:,:,1).^2 + gradSpeed(:,:,:,2).^2 + gradSpeed(:,:,:,3).^2).^0.5;
    parallelity = dot(normal,gradSpeed,4) ./ (normNorm .* normSpeed);
    parallelity(parallelity<parallelityCutoff) = 0;
    
    %% PART 5.2: ADJUST IMAGE INTENSITY
    if parallelitySwitch > 0 && onlyVel < 0
        allEdge = double(allEdge .* parallelity);
    end
    allEdge = int16(allEdge);
    size(allEdge)
    %% PART 5.3: SAVE ALLEDGE TO FILE
    % Edge intensity
    
    %clear normal; clear gradVelocity; clear velocity; clear intensity; clear gradIntensity;
    for z = startZ:1:dimZ
        currSlice = allEdge(:,:,z);
        filename = sprintf('%s%03d%s%03d%s','Cph_',t-1,'_Sec_',z,'.mag');
        filename = strcat(strcat(outDir,'/mag/'), filename);
        %filename = strcat(pwd,filename);
        information{z,1}.StudyID = expt;
        % Note that the information array is inverted
        
        % Try/catch
        worked = -1;
        while worked < 0
            try
                dicomwrite(currSlice,filename,information{z,1},'CreateMode','copy','WritePrivate',true);
                worked = 1;
            catch
                0;
            end
        end
    end
    
    test='Edge-saving complete';
    %clear allEdge;
    %% PART 6 CALCULATE GRAD VELOCITY
    
    % partial matrix of the form shown in the WSS_workflow
    %d/dx
    for ddx=1:1:3
        gradVelocityEdge(:,:,:,ddx,1) = convn(velocity(:,:,:,ddx),xKernelFwd,'same') / delta(1);
        gradVelocityEdge(:,:,:,ddx,2) = convn(velocity(:,:,:,ddx),xKernelRvr,'same') / delta(1);
    end
    % d/dy
    for ddy=1:1:3
        gradVelocityEdge(:,:,:,ddy,3) = convn(velocity(:,:,:,ddy),yKernelFwd,'same') / delta(2);
        gradVelocityEdge(:,:,:,ddy,4) = convn(velocity(:,:,:,ddy),yKernelRvr,'same') / delta(2);
    end
    % d/dz
    for ddz = 1:1:3
        gradVelocityEdge(:,:,:,ddz,5) = convn(velocity(:,:,:,ddz),zKernelFwd,'same') / delta(3);
        gradVelocityEdge(:,:,:,ddz,6) = convn(velocity(:,:,:,ddz),zKernelRvr,'same') / delta(3);
    end
    
    %% PART 6.5 SELECT FOR CORRECT DIRECTION
     for ii = 1:1:size(gradIntensity,1)
         for jj = 1:1:size(gradIntensity,2)
             for kk = 1:1:size(gradIntensity,3)
                 if gradIntensity(ii,jj,kk,1) * delta(1) > 0
                     gradVelocityEdge(ii,jj,kk,:,2) = 0;
                 else
                     gradVelocityEdge(ii,jj,kk,:,1) = 0;
                 end
                 if gradIntensity(ii,jj,kk,2)  * delta(2)> 0
                     gradVelocityEdge(ii,jj,kk,:,4) = 0;
                 else
                     gradVelocityEdge(ii,jj,kk,:,3) = 0;
                 end
                 if gradIntensity(ii,jj,kk,3) * delta(3) > 0
                     gradVelocityEdge(ii,jj,kk,:,6) = 0;
                 else
                     gradVelocityEdge(ii,jj,kk,:,5) = 0;
                 end
             end
         end
     end
    gradVelocity(:,:,:,:,1) = gradVelocityEdge(:,:,:,:,1) + gradVelocityEdge(:,:,:,:,2);
    gradVelocity(:,:,:,:,2) = gradVelocityEdge(:,:,:,:,3) + gradVelocityEdge(:,:,:,:,4);
    gradVelocity(:,:,:,:,3) = gradVelocityEdge(:,:,:,:,5) + gradVelocityEdge(:,:,:,:,6);
    %% PART 7: CALCULATE DIRECTIONAL DERIVATIVE
    % Matrix product of gradIntensity * column vec of normal
    xShear = double(viscosity * (alpha .* gradVelocity(:,:,:,1,1) + beta .* gradVelocity(:,:,:,1,2) + gamma .* gradVelocity(:,:,:,1,3)));
    yShear = double(viscosity * (alpha .* gradVelocity(:,:,:,2,1) + beta .* gradVelocity(:,:,:,2,2) + gamma .* gradVelocity(:,:,:,2,3)));
    zShear = double(viscosity * (alpha .* gradVelocity(:,:,:,3,1) + beta .* gradVelocity(:,:,:,3,2) + gamma .* gradVelocity(:,:,:,3,3)));
    %clear alpha; clear beta; clear gamma; clear gradVelocity;
    if parallelitySwitch > 0 && onlyVel > 0
        xShear = int16(xShear .* parallelity);
        yShear = int16(yShear .* parallelity);
        zShear = int16(zShear .* parallelity);
    else
        xShear = int16(xShear);
        yShear = int16(yShear);
        zShear = int16(zShear);
    end
    
    %% Vel0
    text='saving vel0...';
    
    for z = startZ:1:dimZ
        currSlice = xShear(:,:,z);
        filename = sprintf('%s%03d%s%03d%s','vel_0_corr-',z-1,'-',t-1,'.dcm');
        filename = strcat(strcat(outDir,'/vel_0_corr/'), filename);
        %filename = strcat(pwd,filename);
        % Note that the information array is inverted
        information{z,2}.StudyID = expt;
        % Try/catch
        worked = -1;
        while worked < 0
            try
                dicomwrite(currSlice,filename,information{z,2},'CreateMode','copy','WritePrivate',true);
                worked = 1;
            catch
                0;
            end
        end
    end
    
    %% Vel1
    
    text='saving vel1...';
    
    for z = startZ:1:dimZ
        currSlice = yShear(:,:,z);
        filename = sprintf('%s%03d%s%03d%s','vel_1_corr-',z-1,'-',t-1,'.dcm');
        filename = strcat(strcat(outDir,'/vel_1_corr/'), filename);
        %filename = strcat(pwd,filename);
        % Note that the information array is inverted
        information{z,3}.StudyID = expt;
        % Try/catch
        worked = -1;
        while worked < 0
            try
                dicomwrite(currSlice,filename,information{z,3},'CreateMode','copy','WritePrivate',true);
                worked = 1;
            catch
                0;
            end
        end
    end
    
    %% Vel2
    
    text='saving vel2...';
    
    for z = startZ:1:dimZ
        currSlice = zShear(:,:,z);
        filename = sprintf('%s%03d%s%03d%s','vel_2_corr-',z-1,'-',t-1,'.dcm');
        filename = strcat(strcat(outDir,'/vel_2_corr/'), filename);
        %filename = strcat(pwd,filename);
        % Note that the information array is inverted
        information{z,4}.StudyID = expt;
        % Try/catch
        worked = -1;
        while worked < 0
            try
                dicomwrite(currSlice,filename,information{z,4},'CreateMode','copy','WritePrivate',true);
                worked = 1;
            catch
                0;
            end
        end
    end
    
    text='vel saving done!';
end

end