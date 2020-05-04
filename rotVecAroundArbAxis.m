function rotatedUnitVector = rotVecAroundArbAxis(unitVec2Rotate,rotationAxisUnitVec,theta)
%% Purpose:
%
%  This routine will allow a unit vector, unitVec2Rotate, to be rotated 
%  around an axis defined by the RotationAxisUnitVec.  This is performed by
%  first rotating the unit vector around it's own cartesian axis (in this
%  case we will rotate the vector around the z-axis, [0 0 1]) corresponding
%  to each rotation angle specified by the user via the variable theta ...
%  this rotated vector is then transformed around the user defined axis of
%  rotation as defined by the rotationAxisUnitVec variable.
%  
%
%% References:
%  Murray, G. Rotation About an Arbitrary Axis in 3 Dimensions. 06/06/2013.
%  http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/
%
%% Inputs:
%  unitVec2Rotate               [N x 3]                 Unit Vector in
%                                                       Cartesian 
%                                                       Coordinates to 
%                                                       rotate
%                                                       [x,y,z]
%
%
%  rotationAxisUnitVec          [N x 3]                 Unit Vector with
%                                                       respect to the same 
%                                                       cartesian coordinates 
%                                                       used for unitVec2Rotate
%                                                       [x,y,z]
%
%  theta                        [N x 1]                 Angle in degrees
%                                                       in which to rotate
%                                                       the unitVec2Rotate
%                                                       about the Z-axis
%                                                       before transforming
%                                                       it to the
%                                                       RotateionAxisUnitVec
%                                                       This rotation is
%                                                       counter clockwise
%                                                       when theta is
%                                                       positive, clockwise
%                                                       when theta is
%                                                       negative.
%
%% Outputs:
%  rotatedUnitVector            [N x 3]                 Resulting vector
%                                                       of rotating the
%                                                       unitVec2Rotate
%                                                       about the z-axis
%                                                       described by the
%                                                       angle theta, then
%                                                       transforming the
%                                                       rotated vectors
%                                                       with respect to the
%                                                       rotateionAxisUnitVec
%
%% Revision History:    
%  Darin C. Koblick                                        (c)  03-03-2015
%
%  Darin C. Koblick      Fixed order of rotations               07-30-2015                                       
%% ---------------------- Begin Code Sequence -----------------------------
if nargin == 0
         unitVec2Rotate = [1 0 1]./norm([1 0 1]);
    rotationAxisUnitVec = [1 1 1]./norm([1 1 1]);
    theta = (0:5:360)';
    rotatedUnitVector = rotVecAroundArbAxis(unitVec2Rotate,rotationAxisUnitVec,theta);
    %Show a graphical representation of the rotated vector:
    figure('color',[1 1 1]);
    quiver3(zeros(numel(theta),1),zeros(numel(theta),1),zeros(numel(theta),1), ...
            rotatedUnitVector(:,1),rotatedUnitVector(:,2),rotatedUnitVector(:,3),'k','linewidth',2);
    hold on;
    quiver3(0,0,0,rotationAxisUnitVec(1), ...
                  rotationAxisUnitVec(2), ...
                  rotationAxisUnitVec(3), ...
                  'r','linewidth',5);
    axis equal;
    return;
end
%Check the dimensions of the input vectors to see if we need to repmat
%them:
if size(unitVec2Rotate,1) == 1
   unitVec2Rotate = repmat(unitVec2Rotate,[numel(theta) 1]); 
end
if size(rotationAxisUnitVec,1) == 1
   rotationAxisUnitVec = repmat(rotationAxisUnitVec,[numel(theta) 1]); 
end
%% Step One: take the unit vector rotation axis and rotate into z:
R2Z = vecRotMat(rotationAxisUnitVec,repmat([0 0 1],[size(rotationAxisUnitVec,1) 1]));
unitVectortoRotateAboutZ =Dim33Multiply(unitVec2Rotate,R2Z);
% Rotate the unit vector about the z-axis:
rotatedAboutZAxisUnitVec = bsxRz(unitVectortoRotateAboutZ,theta.*pi/180); 
%% Step Two: Find the rotation Matrix to transform the z-axis to rotationAxisUnitVec
R = vecRotMat(repmat([0 0 1],[size(rotationAxisUnitVec,1) 1]),rotationAxisUnitVec);
%% Step Three: Apply the Rotation matrix to the rotatedAboutZAxisUnitVec vectors
rotatedUnitVector =Dim33Multiply(rotatedAboutZAxisUnitVec,R);
end

function a = Dim33Multiply(a,b)
%% Purpose:
% Given a, an [N x 3] matrix, use b, an [3 x 3 x N] rotation matrix to come
% up with a vectorized solution to b*a
%
%% Inputs:
%  a            [N x 3]                                        N x 3 vector
%
%  b            [3 x 3 x N]                                    3 x 3 x N
%                                                              matrix
%
%% Outputs:
%  a            [N x 3]                                        vectorized
%                                                              solution
%                                                              a = b*a
%
%% Revision History:
% Created by Darin C. Koblick   (C)                                    2013
%% ---------------------- Begin Code Sequence -----------------------------
a =cat(1,sum(permute(bsxfun(@times,b(1,:,:),permute(a,[3 2 1])),[2 3 1]),1), ...
         sum(permute(bsxfun(@times,b(2,:,:),permute(a,[3 2 1])),[2 3 1]),1), ...
         sum(permute(bsxfun(@times,b(3,:,:),permute(a,[3 2 1])),[2 3 1]),1))';
end

function R = vecRotMat(f,t)
%% Purpose:
%Commonly, it is desired to have a rotation matrix which will rotate one 
%unit vector, f,  into another unit vector, t. It is desired to 
%find R(f,t) such that R(f,t)*f = t.  
%
%This program, vecRotMat is the most
%efficient way to accomplish this task. It uses no square roots or
%trigonometric functions as they are very computationally expensive. 
%It is derived from the work performed by Moller and Hughes, which have
%suggested that this method is the faster than any previous transformation
%matrix methods tested.
%
%
%% Inputs:
%f                      [N x 3]                         N number of vectors
%                                                       in which to
%                                                       transform into
%                                                       vector t.
%
%t                      [N x 3]                         N number of vectors
%                                                       in which it is
%                                                       desired to rotate
%                                                       f.
%
%
%% Outputs:
%R                      [3 x 3 x N]                     N number of
%                                                       rotation matrices
%
%% Source:
% Moller,T. Hughes, F. "Efficiently Building a Matrix to Rotate One 
% Vector to Another", 1999. http://www.acm.org/jgt/papers/MollerHughes99
%
%% Created By:
% Darin C. Koblick (C) 07/17/2012
% Darin C. Koblick     04/22/2014       Updated when lines are close to
%                                       parallel by checking 
%% ---------------------- Begin Code Sequence -----------------------------
%It is assumed that both inputs are in vector format N x 3
dim3 = 2;
%Declare function handles for multi-dim operations
normMD = @(x,y) sqrt(sum(x.^2,y));
anyMD  = @(x) any(x(:));
% Inputs Need to be in Unit Vector Format
if anyMD(single(normMD(f,dim3)) ~= single(1)) || anyMD(single(normMD(t,dim3)) ~= single(1))
   error('Input Vectors Must Be Unit Vectors');
end
%Pre-Allocate the 3-D transformation matrix
R = NaN(3,3,size(f,1));

v = permute(cross(f,t,dim3),[3 2 1]);
c = permute(dot(f,t,dim3),[3 2 1]);
h = (1-c)./dot(v,v,dim3);

idx  = abs(c) > 1-1e-13;
%If f and t are not parallel, use the following computation
if any(~idx)
%For any vector u, the rotation matrix is found from:
R(:,:,~idx) = ...
    [c(:,:,~idx) + h(:,:,~idx).*v(:,1,~idx).^2,h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)-v(:,3,~idx),h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)+v(:,2,~idx); ...
     h(:,:,~idx).*v(:,1,~idx).*v(:,2,~idx)+v(:,3,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,2,~idx).^2,h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)-v(:,1,~idx); ...
     h(:,:,~idx).*v(:,1,~idx).*v(:,3,~idx)-v(:,2,~idx),h(:,:,~idx).*v(:,2,~idx).*v(:,3,~idx)+v(:,1,~idx),c(:,:,~idx)+h(:,:,~idx).*v(:,3,~idx).^2];
end
%If f and t are close to parallel, use the following computation
if any(idx)
     f = permute(f,[3 2 1]);
     t = permute(t,[3 2 1]);
     p = zeros(size(f));
     iidx = abs(f(:,1,:)) <= abs(f(:,2,:)) & abs(f(:,1,:)) < abs(f(:,3,:));
     if any(iidx & idx)
        p(:,1,iidx & idx) = 1;
     end
     iidx = abs(f(:,2,:)) < abs(f(:,1,:)) & abs(f(:,2,:)) <= abs(f(:,3,:));
     if any(iidx & idx)
        p(:,2,iidx & idx) = 1;
     end
     iidx = abs(f(:,3,:)) <= abs(f(:,1,:)) & abs(f(:,3,:)) < abs(f(:,2,:));
     if any(iidx & idx)
        p(:,3,iidx & idx) = 1;
     end
     u = p(:,:,idx)-f(:,:,idx);
     v = p(:,:,idx)-t(:,:,idx);
     rt1 = -2./dot(u,u,dim3);
     rt2 = -2./dot(v,v,dim3);
     rt3 = 4.*dot(u,v,dim3)./(dot(u,u,dim3).*dot(v,v,dim3));
     R11 = 1 + rt1.*u(:,1,:).*u(:,1,:)+rt2.*v(:,1,:).*v(:,1,:)+rt3.*v(:,1,:).*u(:,1,:);
     R12 = rt1.*u(:,1,:).*u(:,2,:)+rt2.*v(:,1,:).*v(:,2,:)+rt3.*v(:,1,:).*u(:,2,:);
     R13 = rt1.*u(:,1,:).*u(:,3,:)+rt2.*v(:,1,:).*v(:,3,:)+rt3.*v(:,1,:).*u(:,3,:);
     R21 = rt1.*u(:,2,:).*u(:,1,:)+rt2.*v(:,2,:).*v(:,1,:)+rt3.*v(:,2,:).*u(:,1,:);
     R22 = 1 + rt1.*u(:,2,:).*u(:,2,:)+rt2.*v(:,2,:).*v(:,2,:)+rt3.*v(:,2,:).*u(:,2,:);
     R23 = rt1.*u(:,2,:).*u(:,3,:)+rt2.*v(:,2,:).*v(:,3,:)+rt3.*v(:,2,:).*u(:,3,:);
     R31 = rt1.*u(:,3,:).*u(:,1,:)+rt2.*v(:,3,:).*v(:,1,:)+rt3.*v(:,3,:).*u(:,1,:);
     R32 = rt1.*u(:,3,:).*u(:,2,:)+rt2.*v(:,3,:).*v(:,2,:)+rt3.*v(:,3,:).*u(:,2,:);
     R33 = 1 + rt1.*u(:,3,:).*u(:,3,:)+rt2.*v(:,3,:).*v(:,3,:)+rt3.*v(:,3,:).*u(:,3,:);
     R(:,:,idx) = [R11 R12 R13; R21 R22 R23; R31 R32 R33];
end
end

function m = bsxRz(m,theta)
%% Purpose:
% Perform a rotation of theta radians about the z-axis on the vector(s) 
% described by m.
% 
%% Inputs:
% m         [N x 3]                                         vector matrix
%                                                           in which you
%                                                           would like to
%                                                           rotate with the
%                                                           x,y,z
%                                                           components
%                                                           specified along
%                                                           a specific
%                                                           dimension
%
% theta     [N x 1]                                         Rotation Angle
%                                                           about z-axis 
%                                                           in radians
%
%% Outputs:
% m         [N x 3]
%
%% Revision History:
%  Darin C Koblick (C)                              Initially Created 2013
%% ---------------------- Begin Code Sequence -----------
%Assemble the rotation matrix
Rz = zeros(3,3,size(m,1));
Rz(1,1,:) = cos(theta);  Rz(1,2,:) = -sin(theta);
Rz(2,1,:) = sin(theta);  Rz(2,2,:) =  cos(theta); 
Rz(3,3,:) = 1;
%Dim33Multiply
m = Dim33Multiply(m,Rz);
end