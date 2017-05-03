1;
% newamp1_getvqset.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License
% Version 1

pkg load statistics

function surface_vector = generate_save_training_set(input_prefix, output_prefix)
    newamp;
    more off;

    max_amp = 80;
    postfilter = 0;   % optional postfiler that runs on Am, not used atm
    synth_phase = 1;

    #if nargin == 1
    #  output_prefix = input_prefix;
    #end
    model_name = strcat(input_prefix,"_model.txt");
    model = load(model_name);
    [frames nc] = size(model);
    surface_no_mean = rate_K_dec_vq_dump_xfbf(model);
    surface_name = strcat(output_prefix,"_vectors.h5");
    save("-hdf5",surface_name,"surface_no_mean")
    surface_vector = surface_no_mean;
endfunction


function surface = rate_K_dec_vq_dump_xfbf(model)
    M = 4; % model frame -> wire frame decimation rate 10ms->40ms
    K = 20;

    max_amp = 80;
    [frames nc] = size(model);
    xframes = floor(frames/M)
    surface = zeros(xframes,K);
    Wof = zeros(xframes);
    mean_f_i = zeros(xframes);

    %Surface prediction coeffecient
    coeff_surface_pred = .8;
    frame_last = zeros(1,K);
    melvq;
    % create frames x K surface.  TODO make all of this operate frame by
    % frame, or at least M/2=4 frames rather than one big chunk

    %Do the compression part, frame by frame
    for fx = 1:xframes
        f = ((fx-1)*M)+1;
        model_frame = model(f,:);
        [frame_mel sample_kHz] = resample_const_rate_f_mel(model_frame,K);
        %Calculate and remove mean (in dB)
        mean_f = mean(frame_mel);
        frame_no_mean = frame_mel - mean_f;
        %Remove overall slope
        b = regress(frame_no_mean',sample_kHz');
        n = sample_kHz*b;

        %Not removing or regenerating mean for this test
        frame_no_slope = frame_mel ;
        %Figure out what the other side expects this frame to look like
        frame_pred = frame_last .* coeff_surface_pred;
        frame_resid = frame_no_slope - frame_pred;

        %[res frame_resid_ ind] = mbest(train_120_vq, frame_resid, m);

        frame_resid_ = frame_resid;
        surface(fx,:) = frame_resid_;
        frame_last = frame_pred + frame_resid;

        printf("\rEncode frame %d of %d (%d/%d)",f,frames,fx,xframes);
    end
    printf("\n");

endfunction

 % -----------------------------------------------------------------------------------------
 % Linear decimator/interpolator that operates at rate K, includes VQ, post filter, and Wo/E
 % quantisation.  Evolved from abys decimator below.  Simulates the entire encoder/decoder.

% -----------------------------------------------------------------------------------------
% Linear decimator/interpolator that operates at rate K, includes VQ, post filter, and Wo/E
% quantisation.  Evolved from abys decimator below.  Simulates the entire encoder/decoder.
%
% Stripped down to dump the rate-k surface as a VQ training set

function surface_no_mean = rate_K_dec_vq_dump(model)
  max_amp = 80;
  [frames nc] = size(model);
  model_ = zeros(frames, max_amp+3);
  indexes = zeros(frames,4);

  M = 4;

  % create frames x K surface.  TODO make all of this operate frame by
  % frame, or at least M/2=4 frames rather than one big chunk

  K = 20;
  [surface sample_freqs_kHz] = resample_const_rate_f_mel(model, K);
  target_surface = surface;

  %figure(1);
  %mesh(surface);

  % VQ rate K surface.  TODO: If we are decimating by M/2=4 we really
  % only need to do this every 4th frame.

  melvq;
  load train_120_vq; m=5;

  for f=1:frames
    mean_f(f) = mean(surface(f,:));
    surface_no_mean(f,:) = surface(f,:) - mean_f(f);
    b = regress(surface_no_mean(f,:)',sample_freqs_kHz');
    n = sample_freqs_kHz*b;
    surface_no_mean(f,:) = surface(f,:) - n - mean_f(f);
    printf("\rrm mean frame %d of %d",f,frames)
  end
  printf("\n")
  figure(1)
  hist(mean_f)
  %figure(2)

  %figure(2);
  %mesh(surface_no_mean);
endfunction
