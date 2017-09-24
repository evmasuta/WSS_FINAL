function output = KernelMakerv4(diff,space,sigma)
% Returns, in order, the kernels for d/dx, d/dy, d/dz
output = cell(4);
% unpack spacings
dr = space(1);
dc = space(2);
ds = space(3);
% Kernel1 == d/dr; Kernel2 == d/dc; Kernel3 == d/ds
%% Generate kernel1
kernel1 = zeros([5 5 5]);
kernel4 = kernel1;
kernel1(:,3,3) = fliplr([0,0,-3,4,-1]) ./ 2;
kernel4(:,3,3) = fliplr([1,-4,3,0,0]) ./ 2;
%% Generate kernel2
kernel2 = zeros([5 5 5]);
kernel5 = kernel2;
kernel2(3,:,3) = fliplr([0,0,-3,4,-1]) ./ 2;
kernel5(3,:,3) = fliplr([1,-4,3,0,0]) ./ 2;
%% Generate kernel3
kernel3 = zeros([5 5 5]);
kernel6 = kernel3;
kernel3(3,3,:) = fliplr([0,0,-3,4,-1]) ./ 2;
kernel6(3,3,:) = fliplr([1,-4,3,0,0]) ./ 2;
tempCell = {kernel1, kernel2, kernel3};
tempCell2 = {kernel4,kernel5,kernel6};
% Determine indices
xIndex = find(diff(:,1));
yIndex = find(diff(:,2));
zIndex = find(diff(:,3));
spacings(1) = diff(xIndex,1);
spacings(2) = diff(yIndex,2);
spacings(3) = diff(zIndex,3);
% NEED TO DIVIDE BY 10 TO CONVERT FROM mm to cm!!!
spacings = spacings / 10;
% Sort kernels to xKernel, yKernel, zKernel, -xKernel, -yKernel, -zKernel

output{1} = tempCell{xIndex};
output{2} = tempCell{yIndex};
output{3} = tempCell{zIndex};
output{4} = tempCell2{xIndex};
output{5} = tempCell2{yIndex};
output{6} = tempCell2{zIndex};
output{7} = spacings;