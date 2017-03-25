%

fsk_horus_as_a_lib = 1; % make sure calls to test functions at bottom are disabled

fsk_horus
pkg load signal;
pkg load parallel;
graphics_toolkit('gnuplot');


function test_stats = fsk_demod_xt(Fs,Rs,f1,fsp,mod,tname,M=2)
    o_f1_dc = [];
    o_f2_dc = [];
    o_f3_dc = [];
    o_f4_dc = [];
    o_f1_int = [];
    o_f2_int = [];
    o_f3_int = [];
    o_f4_int = [];
    o_f1 = [];
    o_f2 = [];
    o_f3 = [];
    o_f4 = [];
    o_EbNodB = [];
    o_ppm = [];
    o_Sf = [];
    o_fest = [];
    o_rx_timing = [];
    o_norm_rx_timing = [];
    o_nin = [];
    %Run octave demod, dump some test vectors
    states = fsk_horus_init_hbr(Fs,10,Rs,M);
    Ts = states.Ts;
    P = states.P;
    states.ftx(1) = f1;
    states.ftx(2) = f1+fsp;
    states.ftx(3) = f1+fsp*2;
    states.ftx(4) = f1+fsp*3;
    states.dA = 1;
    states.dF = 0;
    modin = mod;
    obits = [];
    while length(modin)>=states.nin
        ninold = states.nin;
        states = est_freq(states, modin(1:states.nin), states.M);
        [bitbuf,states] = fsk_horus_demod(states, modin(1:states.nin));
        modin=modin(ninold+1:length(modin));
        obits = [obits bitbuf];

        %Save other parameters
        o_f1_dc = [o_f1_dc states.f_dc(1,1:states.Nmem-Ts/P)];
        o_f2_dc = [o_f2_dc states.f_dc(2,1:states.Nmem-Ts/P)];
        o_f1_int = [o_f1_int states.f_int(1,:)];
        o_f2_int = [o_f2_int states.f_int(2,:)];
        o_EbNodB = [o_EbNodB states.EbNodB];
        o_ppm = [o_ppm states.ppm];
        o_rx_timing = [o_rx_timing states.rx_timing];
        o_norm_rx_timing = [o_norm_rx_timing states.norm_rx_timing];
        o_Sf = [o_Sf states.Sf'];
        o_f1 = [o_f1 states.f(1)];
        o_f2 = [o_f1 states.f(2)];
        o_fest = [o_fest states.f];
        o_nin = [o_nin states.nin];
        if M==4
      			o_f3_dc = [o_f3_dc states.f_dc(3,1:states.Nmem-Ts/P)];
      			o_f4_dc = [o_f4_dc states.f_dc(4,1:states.Nmem-Ts/P)];
      			o_f3_int = [o_f3_int states.f_int(3,:)];
      			o_f4_int = [o_f4_int states.f_int(4,:)];
      			o_f3 = [o_f1 states.f(3)];
      			o_f4 = [o_f1 states.f(4)];
        end
    end

    %close all

    pass = 1;

    test_stats.pass = pass;
    test_stats.obits = obits;

endfunction


% A big ol' channel impairment tester
% Shamlessly taken from fsk_horus
% This throws some channel imparment or another at the C and octave modem so they
% may be compared.
function stats = tfsk_run_sim(test_frame_mode,EbNodB,timing_offset,fading,df,dA,M=2)
  global print_verbose;
  frames = 250;
  %EbNodB = 10;
  %timing_offset = 2.0; % see resample() for clock offset below
  %fading = 0;          % modulates tx power at 2Hz with 20dB fade depth,
                       % to simulate balloon rotating at end of mission
  %df     = 0;          % tx tone freq drift in Hz/s
  %dA     = 1;          % amplitude imbalance of tones (note this affects Eb so not a gd idea)

  more off
  rand('state',1);
  randn('state',1);

  % ----------------------------------------------------------------------

  % sm2000 config ------------------------
  %states = fsk_horus_init(96000, 1200);
  %states.f1_tx = 4000;
  %states.f2_tx = 5200;

  if test_frame_mode == 2
    % 2400A config
    states = fsk_horus_init_hbr(48000,10, 1200, M);
    states.f1_tx = 1200;
    states.f2_tx = 2400;
    states.f3_tx = 3600;
    states.f4_tx = 4800;
    states.ftx(1) = 1200;
    states.ftx(2) = 2400;
    states.ftx(3) = 3600;
    states.ftx(4) = 4800;

  end

  if test_frame_mode == 4
    % horus rtty config ---------------------
    states = fsk_horus_init_hbr(48000,10, 1200, M);
    states.f1_tx = 1200;
    states.f2_tx = 2400;
    states.f3_tx = 3600;
    states.f4_tx = 4800;
    states.ftx(1) = 1200;
    states.ftx(2) = 2400;
    states.ftx(3) = 3600;
    states.ftx(4) = 4800;

    states.tx_bits_file = "horus_tx_bits_rtty.txt"; % Octave file of bits we FSK modulate

  end

  if test_frame_mode == 5
    % horus binary config ---------------------
    states = fsk_horus_init_hbr(48000,10, 1200, M);
    states.f1_tx = 1200;
    states.f2_tx = 2400;
    states.f3_tx = 3600;
    states.f4_tx = 4800;
    states.ftx(1) = 1200;
    states.ftx(2) = 2400;
    states.ftx(3) = 3600;
    states.ftx(4) = 4800;
    %%%states.tx_bits_file = "horus_tx_bits_binary.txt"; % Octave file of bits we FSK modulate
	states.tx_bits_file = "horus_payload_rtty.txt";
  end

  % ----------------------------------------------------------------------

  states.verbose = 0;%x1;
  N = states.N;
  P = states.P;
  Rs = states.Rs;
  nsym = states.nsym;
  Fs = states.Fs;
  states.df = df;
  states.dA = [dA dA dA dA];
  states.M = M;

  EbNo = 10^(EbNodB/10);
  variance = states.Fs/(states.Rs*EbNo*states.bitspersymbol);

  % set up tx signal with payload bits based on test mode

  if test_frame_mode == 1
     % test frame of bits, which we repeat for convenience when BER testing
    test_frame = round(rand(1, states.nsym));
    tx_bits = [];
    for i=1:frames+1
      tx_bits = [tx_bits test_frame];
    end
  end
  if test_frame_mode == 2
    % random bits, just to make sure sync algs work on random data
    tx_bits = round(rand(1, states.nbit*(frames+1)));
  end
  if test_frame_mode == 3
    % ...10101... sequence
    tx_bits = zeros(1, states.nsym*(frames+1));
    tx_bits(1:2:length(tx_bits)) = 1;
  end

  if (test_frame_mode == 4) || (test_frame_mode == 5)

    % load up a horus msg from disk and modulate that

    test_frame = load(states.tx_bits_file);
    ltf = length(test_frame);
    ntest_frames = ceil((frames+1)*nsym/ltf);
    tx_bits = [];
    for i=1:ntest_frames
      tx_bits = [tx_bits test_frame];
    end
  end



  f1 = states.f1_tx;
  fsp = states.f2_tx-f1;
  states.dA = [dA dA dA dA];
  states.ftx(1) = f1;
  states.ftx(2) = f1+fsp;

  if states.M == 4
	states.ftx(3) = f1+fsp*2;
	states.ftx(4) = f1+fsp*3;
  end

  tx = fsk_horus_mod(states, tx_bits);

  if timing_offset
    tx = resample(tx, 1000, 1001); % simulated 1000ppm sample clock offset
  end

  if fading
     ltx = length(tx);
     tx = tx .* (1.1 + cos(2*pi*2*(0:ltx-1)/Fs))'; % min amplitude 0.1, -20dB fade, max 3dB
  end

  noise = sqrt(variance)*randn(length(tx),1);

  %chirp_f0
  %chirp_rate =

  % One full chirp per symbol
  Ts = Fs/Rs;
  chirp_f0 = 2*pi*(0/Fs);
  chirp_f1 = 2*pi*(4800/Fs);
  chirp_k = (chirp_f1-chirp_f0)/(Ts);
  %chirp_t = (1:Ts);
  %chirp_single = sin(2*pi*( chirp_f0*chirp_t + (chirp_k/2)*chirp_t.^2 )/Fs);

  chirp_th = 0; % Chirp Theta
  chirp_fl = (1:Ts)*chirp_k;


  chirp_out = zeros(1,length(tx));
  for ii = 1:length(tx)
    chirp_f_i = mod(ii,Ts);
    chirp_f = chirp_fl(chirp_f_i+1);
    chirp_out(ii) = exp(j*chirp_th);
    chirp_th = chirp_th + chirp_f;
  end
  chirp_out = chirp_out';
  %plot(chirp_out)
  tx_chirp = tx .* chirp_out;
  %tx_chrip = real(tx_chirp);

  rx_chirp = tx_chirp + noise;

  rx = rx_chirp .* conj(chirp_out);
  rx = real(rx);
  %plot(20*log10(abs(fft(rx))))

%rx    = tx + noise;

  test_name = sprintf("tfsk EbNodB:%d frames:%d timing_offset:%d fading:%d df:%d",EbNodB,frames,timing_offset,fading,df);
  tstats = fsk_demod_xt(Fs,Rs,states.f1_tx,fsp,rx,test_name,M);

  pass = tstats.pass;
  obits = tstats.obits;
  %cbits = tstats.cbits;
  %stats.name = test_name;

  %if tstats.pass
  %  printf("Test %s passed\n",test_name);
  %else
  %  printf("Test %s failed\n",test_name);
  %end

  % Figure out BER of octave and C modems
  bitcnt = length(tx_bits);
  rx_bits = obits;
  ber = 1;
  ox = 1;
  for offset = (1:100)
    nerr = sum(xor(rx_bits(offset:length(rx_bits)),tx_bits(1:length(rx_bits)+1-offset)));
    bern = nerr/(bitcnt-offset);
    if(bern < ber)
      ox = offset;
      best_nerr = nerr;
    end
    ber = min([ber bern]);
  end
  offset = ox;
  bero = ber;

  offset = ox;
  stats.ber = bero;

  % coherent BER theory calculation
  stats.thrcoh = .5*(M-1)*erfc(sqrt( (log2(M)/2) * EbNo ));

  % non-coherent BER theory calculation
  % It was complicated, so I broke it up
  ms = M;
  ns = (1:ms-1);
  as = (-1).^(ns+1);
  bs = (as./(ns+1));

  cs = ((ms-1)./ns);

  ds = ns.*log2(ms);
  es = ns+1;
  fs = exp( -(ds./es)*EbNo );

  thrncoh = ((ms/2)/(ms-1)) * sum(bs.*((ms-1)./ns).*exp( -(ds./es)*EbNo ));

  stats.thrncoh = thrncoh;
  stats.pass = pass;
endfunction


function plot_fsk_bers(M=2)
    %Range of EbNodB over which to plot
    ebnodbrange = (4:13);

    %berc = ones(1,length(ebnodbrange));
    bero = ones(1,length(ebnodbrange));
    berinc = ones(1,length(ebnodbrange));
    beric = ones(1,length(ebnodbrange));
    ebnodbs = length(ebnodbrange)
    mode = 2;
    %Replication of other parameters for parcellfun
    modev   = repmat(mode,1,ebnodbs);
    timingv = repmat(1,1,ebnodbs);
    fadingv = repmat(0,1,ebnodbs);
    dfv     = repmat(1,1,ebnodbs);
    dav     = repmat(1,1,ebnodbs);
    Mv     = repmat(M,1,ebnodbs);


    statv = pararrayfun(floor(nproc()),@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav,Mv);
    %statv = arrayfun(@tfsk_run_sim,modev,ebnodbrange,timingv,fadingv,dfv,dav,Mv);

    for ii = (1:length(statv))
        stat = statv(ii);
        %berc(ii)=stat.berc;
        bero(ii)=stat.ber;
        berinc(ii)=stat.thrncoh;
        beric(ii) = stat.thrcoh;
    end
    clf;
    figure(M)

    semilogy(ebnodbrange, berinc,sprintf('r;%dFSK non-coherent theory;',M))
    hold on;
    semilogy(ebnodbrange, beric ,sprintf('g;%dFSK coherent theory;',M))
    semilogy(ebnodbrange, bero ,sprintf('b;Octave fsk horus %dFSK Demod;',M))
    %semilogy(ebnodbrange, berc,sprintf('+;C fsk horus %dFSK Demod;',M))
    hold off;
    grid("minor");
    axis([min(ebnodbrange) max(ebnodbrange) 1E-5 1])
    legend("boxoff");
    xlabel("Eb/No (dB)");
    ylabel("Bit Error Rate (BER)")

endfunction
