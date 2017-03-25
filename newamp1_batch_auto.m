% newamp1_batch.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License
% Version 1

%This program is a helper to automate codec idea testing with newamp1_batch
% It automatically runs c2sim and newamp1_batch to go from raw->model->raw


function newamp1_batch_auto(raw_file_prefix_in,raw_file_prefix_out):
    %Path to c2sim binary. Change to suit your system
    c2sim_path = "../build_linux/src/c2sim"
    pre = raw_file_prefix_in
    command = sprintf("%s %s.raw --dump %s --phase0 --lpc 10 --dump_pitch_e %s",c2sim_path,pre,pre,pre)
    system(command)
    newamp1_batch(pre,raw_file_prefix_out)
    ipre = raw_file_prefix_in
    pre = raw_file_prefix_out
    command = sprintf("%s %s.raw --phase0 --postfilter --amread %s_am.out --hmread %s_hm.out --Woread %s_Wo.out -o %s_out.raw",c2sim_path,ipre,pre,pre,pre,pre)
