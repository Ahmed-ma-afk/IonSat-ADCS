% Function to create the point cloud of a basic CubeSat. 
%
% -----------------------------------------------------------------------------
%   Copyright (c) 2010-2018 Samir A. Rawashdeh
%   Electrical and Computer Engineering
%   University of Michigan - Dearborn
%  
%   All rights reserved. 
%   
%   Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are 
%   met:
%   
%       * Redistributions of source code must retain the above copyright 
%         notice, this list of conditions and the following disclaimer.
%       * Redistributions in binary form must reproduce the above copyright 
%         notice, this list of conditions and the following disclaimer in 
%         the documentation and/or other materials provided with the distribution
%         
%   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%   POSSIBILITY OF SUCH DAMAGE.
%  
% ----------------------------------------------------------------------------

function [cg, volume] = draw_cubesatv1(body_length, cg_offset, volume, Resolution)

volume = volume * NaN;
len = body_length / Resolution;  % length in dots = 34cm * Resolution
%Original 3U
%wid = 10 / Resolution;  % width in dots = 10cm * Resolution
%Modified MDQsat-1
%wid = 5 / Resolution;  % width in dots = 10cm * Resolution
%Modified IonSat
wid = 22 / Resolution;  % width in dots = 22.63cm * Resolution (can only be integer), before 23, now 22
hei = 10 / Resolution;  % width in dots = 10cm * Resolution
cg = round(size(volume)/2); % temporary cg - center of the volume
%cg = [180 180 180];



%% main satellite body
x_range = (cg-len/2):(cg+len/2); %112 to 248
y_range = (cg-wid/2):(cg+wid/2); %136 to 224
z_range = (cg-hei/2):(cg+hei/2); %160 to 200
volume(x_range,y_range,z_range) = 1;

%Deployable Solar Panels 
x_range = 112:248;  %180 +- 17cm(/0.25 which is resolution)
y_range = 72:288;   %180 +- 27cm(/0.25 which is resolution)
z_range = 199;      %180 + 5(/0.25 which is resolution) -1 
volume(x_range,y_range,z_range) = 1;

% %Tail, just to test
% x_range = 249:259;  %180 +- 17cm(/0.25 which is resolution)
% y_range = 170:190;   %180 +- 27cm(/0.25 which is resolution)
% z_range = 199;      %180 + 5(/0.25 which is resolution) -1 
% volume(x_range,y_range,z_range) = 1;


%Here we substract the points that are not a surface
x_range = (cg-len/2)+1:(cg+len/2)-1;    %113 to 247
y_range = (cg-wid/2)+1:(cg+wid/2)-1;    %137 to 223
z_range = (cg-hei/2)+1:(cg+hei/2)-1;    %161 to 199
volume(x_range,y_range,z_range) = NaN;

% Example, to add deployables, like antennae or booms:
% x_range = 200:220;
% y_range = 189:190;
% z_range = 189;
% volume(x_range,y_range,z_range) = 1

%Added for IonSat deployable solar panels
%x_range = 148:212; %OK
%y_range = 120:240;
%y_range = 125:235; %OK
%z_range = 188; %OK



%% CG Correction
%% IMPORTANT:  Best to hardcode CG outside of this function, variable called "origin"
% [x y z] = find_1s_in_volume(volume);
% cg = round([sum(x)/length(x), sum(y)/length(y), sum(z)/length(z)]);  % assumes satellite is uniformly distributed
% cg(1) = cg(1) + cg_offset/Resolution;  


        
       



