@@ -5,6 +5,8 @@
 % This program is distributed under the terms of the GNU General Public License
 % Version 1
 
+pkg load statistics
+
 function surface_vector = generate_save_training_set(input_prefix, output_prefix)
     newamp;
     more off;
@@ -19,12 +21,76 @@ function surface_vector = generate_save_training_set(input_prefix, output_prefix
     model_name = strcat(input_prefix,"_model.txt");
     model = load(model_name);
     [frames nc] = size(model);
-    surface_no_mean = rate_K_dec_vq_dump(model);
+    surface_no_mean = rate_K_dec_vq_dump_xfbf(model);
     surface_name = strcat(output_prefix,"_vectors.h5");
     save("-hdf5",surface_name,"surface_no_mean")
     surface_vector = surface_no_mean;
 endfunction
 
+function surface = rate_K_dec_vq_dump_xfbf(model)
+    M = 4; % model frame -> wire frame decimation rate 10ms->40ms
+    MP = M/2; % Multi frame packing rate
+    K = 20;
+
+    max_amp = 80;
+    [frames nc] = size(model);
+    xframes = floor(frames/M)
+    indexes = zeros(xframes,4);
+    surface = zeros(xframes,K*2);
+    Wof = zeros(xframes);
+    mean_f_i = zeros(xframes);
+
+    melvq;
+    % create frames x K surface.  TODO make all of this operate frame by
+    % frame, or at least M/2=4 frames rather than one big chunk
+
+    energy_q = create_energy_q;
+    %Do the compression part, frame by frame
+
+    for fx=1:xframes
+        f = ((fx-1)*M)+1;
+        model_frame_a = model(f,:);
+        model_frame_b = model(f+MP,:);
+        [frame_mel_a sample_kHz] = resample_const_rate_f_mel(model_frame_a,K);
+        [frame_mel_b sample_kHz] = resample_const_rate_f_mel(model_frame_b,K);
+
+        %Calculate and remove mean (in dB)
+        mean_f_a = mean(frame_mel_a);
+        frame_no_mean_a = frame_mel_a - mean_f_a;
+        %Remove overall slope
+        b = regress(frame_no_mean_a',sample_kHz');
+        n = sample_kHz*b;
+        frame_mel2_a = frame_mel_a - n;
+        frame_no_mean_a = frame_mel2_a - mean_f_a;
+
+        %Calculate and remove mean (in dB) from middle frame
+        mean_f_b = mean(frame_mel_b);
+        frame_no_mean_b = frame_mel_b - mean_f_b;
+        %Remove overall slope
+        b = regress(frame_no_mean_b',sample_kHz');
+        n = sample_kHz*b;
+        frame_mel2_b = frame_mel_b - n;
+
+        %Subtract 'A' mean
+        frame_no_mean_b = frame_mel2_b - mean_f_a;
+
+        %Take delta of frame
+        frame_b_dt = frame_mel2_b - frame_mel2_a;
+
+        frame_no_mean = zeros(1,K*2);
+        frame_no_mean(1:2:K*2) = frame_no_mean_a;
+        frame_no_mean(2:2:K*2) = frame_no_mean_b;
+
+        %[res frame_no_mean_ ind] = mbest(train_120_vq, frame_no_mean, m);
+        %indexes(fx,1:2) = ind;
+        frame_no_mean_ = frame_no_mean;
+        surface(fx,:) = frame_no_mean_;
+
+        printf("\rEncode frame %d of %d",f,frames);
+    end
+    printf("\n");
+endfunction
+
 % -----------------------------------------------------------------------------------------
 % Linear decimator/interpolator that operates at rate K, includes VQ, post filter, and Wo/E
 % quantisation.  Evolved from abys decimator below.  Simulates the entire encoder/decoder.
