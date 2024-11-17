% Calculate the CG.
%   CG Calculator for StratoShield, AERE 4610 Fall 2024. All measurements
%   of distance are in inches and weight is in pounds. The position of each
%   part is measured from the tip of the aircraft's nose.
clear; clc; close all;

%% Aircraft values

c_wing = 8; % [in]
c_tail = 6; % [in]
boom_length = 36; % [in]

%% Parts

% Wing foam
W_wing_foam = 0.368; % [lb]
x_wing_foam = 1.5 * c_wing; % [in]

% Wing spar (2 for the wing)
W_wing_spar = 0.2125; % lb
x_wing_spar = 1.5 * c_wing; % in

% Tail foam
W_tail_foam = 0.05; % lb
x_tail_foam = boom_length - c_tail / 2;

% Tail spar (6 for the tail)
W_tail_spar = 0.024; % lb
x_tail_spar = boom_length - c_tail / 2;

% Boom
W_boom = 0.11 / 36 * boom_length; % lb
x_boom = boom_length / 2; % in

% Fuselage
W_fuselage = 1.07625; % lb
x_fuselage = 27 / 2; % in

% Payload
W_payload = 2.8108938; % lb
x_payload = 27 / 2; % in

%% Calculations

W_total = W_wing_foam ...
    + 2 * W_wing_spar ...
    + W_tail_foam ...
    + 6 * W_tail_spar ...
    + W_boom ...
    + W_fuselage ...
    + W_payload;

CG = (W_wing_foam * x_wing_foam ...
    + 2 * W_wing_spar * x_wing_spar ...
    + W_tail_foam * x_tail_foam ...
    + 6 * W_tail_spar * x_tail_spar ...
    + W_boom * x_boom ...
    + W_fuselage * x_fuselage ...
    + W_payload * x_payload) ...
    / W_total;

%% Output

fprintf("Total weight = %g lb\n", W_total);
fprintf("CG position = %g in\n", CG);
