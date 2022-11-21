% Quaternion Library for Simulink
% Version 1.7 (JASP) 12-Dec-2009
% ========================================================================
% Copyright (c) 2000-2009, Jay A. St. Pierre.  All rights reserved.
% This software is licensed under the terms of the BSD license.
% See the file license.txt that came with the software for more details.
% ========================================================================
%
% This is a library of blocks that allows manipulation of quaternions.
%
% The available blocks are:
%
%   Quaternion Normalize
%   Quaternion Conjugate
%   Quaternion Multiply
%
%   Quaternion Propagation
%   Quaternion Vector Transform
%   Quaternion Vector Rotation
%
%   Quaternion Decomposition
%   Quaternion to DCM
%   DCM to Quaternion
%
%   Row Major to Matrix
%   Matrix to Row Major
%
% For purposes of this library, a quaternion, q, is just a four element
% vector where q(1:3) is the "imaginary" or "vector" portion of the
% hypercomplex number, and q(4) is the "real" or "scalar" portion.
% Consequently, if q represents a rotation, then:
%
%   q(1) = v1*sin(phi/2)
%   q(2) = v2*sin(phi/2)
%   q(3) = v3*sin(phi/2)
%   q(4) =    cos(phi/2)
%
% where phi is the amount of rotation about the unit vector [v1 v2 v3].
%
% The DCM's produced by the "Quaternion to DCM" block and used by the DCM
% to Quaternion block are written in row major form instead of normal
% MATLAB matrices.  Back when this library was first developed, Simulink
% signals could not be matrices, so were written in the C-friendly
% row-major format.  To prevent the breakage of models that were written
% to use the row major signals, the interface to the blocks have not been
% changed but two blocks were added to aid with the translation of
% matrices to/from row major format.
%
% IMPORTANT NOTE: For purposes of quaternion/DCM equivalence, the
% relationship is chosen to be:
%
%    R v = q* v q
%
% Therefore a "transform" is (q* v q) and a "rotation" is (q v q*). This
% follows the convention used in "Spacecraft Attitude Determination and
% Control", edited by James R. Wertz, Copyright 1978.  Note that many
% recent uses of quaternions choose the opposite convention (the "left"
% quaternion being the equivalent), including many computer graphics
% libraries.
%
% See also QUATERNIONS, the quaternion manipulation toolbox for the
% MATLAB command line.  RM2MAT and MAT2RM functions provided by the
% Matrix Library for Simulink are also useful for handling
% row-major form matrices.

% Package: $Name: qlib-1_7 $
% File: $Revision: 1.14 $
% $Date: 2009-12-12 20:55:47 $
