1;
% newamp1_getvqset.m
%
% Copyright David Rowe 2017
% This program is distributed under the terms of the GNU General Public License
% Version 1

pkg load statistics
graphics_toolkit('gnuplot')

function [surface wo_e_vec] = generate_save_training_set(input_prefix, output_prefix)
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
    voicing_name = strcat(input_prefix,"_pitche.txt");
    voicing = zeros(1,frames);

    if exist(voicing_name, "file") == 2
      pitche = load(voicing_name);
      voicing = pitche(:, 3);
    end

    [surface wo_e_vec] = experiment_rate_K_dec_xfbf(model,voicing);
    surface_name = strcat(output_prefix,"_vectors.h5");
    save("-hdf5",surface_name,"surface","wo_e_vec")
endfunction

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

function [surface Wo_mean_vec] = experiment_rate_K_dec_xfbf(model,voicing)
  M = 4; % model frame -> wire frame decimation rate 10ms->40ms
  MP = 2; % Multi frame packing rate
  K = 20;

  max_amp = 80;
  [frames nc] = size(model);
  xframes = floor(frames/M); % Frames at rate 1/M
  indexes = zeros(xframes,4);
  surface = zeros(xframes,K);
  Wof = zeros(xframes,1);
  mean_f_i = zeros(xframes,1);
  Wo_mean_vec = zeros(xframes,2);

  melvq;
  % create frames x K surface.  TODO make all of this operate frame by
  % frame, or at least M/2=4 frames rather than one big chunk

  energy_q = create_energy_q;
  last_Wo = 0;
  last_mean = 0;
  mean_predict = .9;
  Wo_predict = .8;
  %Do the compression part, frame by frame
  for fx=1:xframes
      f = ((fx-1)*M)+1;
      model_frame = model(f,:);
      [frame_mel sample_kHz] = resample_const_rate_f_mel(model_frame,K);
      %Calculate and remove mean (in dB)
      mean_f = mean(frame_mel);
      frame_no_mean = frame_mel - mean_f;
      %Remove overall slope
      b = regress(frame_no_mean',sample_kHz');
      n = sample_kHz*b;
      frame_no_mean = frame_mel - n - mean_f;

      %[res frame_no_mean_ ind] = mbest(train_120_vq, frame_no_mean, m);
      %indexes(fx,1:2) = ind;
      frame_no_mean_ = frame_no_mean;
      surface(fx,:) = frame_no_mean_;

      mean_p = last_mean * mean_predict;
      mean_err = mean_p - mean_f;
      Wo_p = last_Wo * Wo_predict;

      %Wo_f = 2*pi/100;
      %if voicing(f)
        Wo_f = model(f,1);
      %end
      Wo_f_log = log10(Wo_f);
      Wo_err = Wo_p - Wo_f_log;

      %Do VQ of mean/Wo here
      Wo_err_ = Wo_err;
      mean_err_ = mean_err;
      Wo_mean_vec(fx,:) = [Wo_err_,mean_err_];

      %Update 'last' values to reflect prediction and VQ error
      last_Wo = Wo_p - Wo_err_;
      last_mean = mean_p - mean_err_;

      %[mean_f_ ind] = quantise(energy_q, mean_f);
      %indexes(fx,3) = ind - 1;
      %mean_f_i(fx) = mean_f_;
      mean_f_i(fx) = mean_f;

      %if voicing(f)
        %ind = encode_log_Wo(model(f,1), 6);
    %    if ind == 0
    %      ind = 1;
    %    end
        %indexes(fx,4) = ind;
        %Wof(fx) = decode_log_Wo(indexes(fx,4), 6);
    %    Wof(fx) = model(f,1);
     % else
    %    indexes(fx,4) = 0;
    %    Wof(fx) = 2*pi/100;
     % end
      printf("\rEncode frame %d of %d",f,frames);
  end
  Wo_mean_vec_t = Wo_mean_vec;
  Wo_mean_vec_t(:,1) = Wo_mean_vec_t(:,1).*80;
  wo_e_vq = trainvq(Wo_mean_vec_t,2**9,1);
  wo_e_vq(:,1) = wo_e_vq(:,1).*(1/80);

  figure(1)
  clf
  hold on
  scatter(Wo_mean_vec(:,1),Wo_mean_vec(:,2),'b','.')
  scatter(wo_e_vq(:,1),wo_e_vq(:,2),10,'r','o','filled')
  ylabel("Energy residual (dB)")
  xlabel("Pitch residual (log10(Wo))")
  hold off
  printf("\n");
end
