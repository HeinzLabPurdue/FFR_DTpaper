function tc_new= get_shifted_tc(tc, bf_new, thresh_new, hearing_status)

% hardcoded for the two units
bf_old= 4;

bf_new_kHz= bf_new;
[~, bf_ind]= min(tc.TCfit);
thresh_old= tc.TCfit(bf_ind);

tc_new= tc;

if strcmp(hearing_status, 'nh')
    tc_new.freqkHz= (tc.freqkHz - bf_old)*bf_new_kHz/bf_old + bf_new_kHz;
elseif strcmp(hearing_status, 'hi')
    tc_new.freqkHz= bf_new_kHz + (tc.freqkHz - bf_old)*(bf_new_kHz - 0.1)/(bf_old - min(tc.freqkHz));
end

tc_new.TCfit= tc.TCfit - thresh_old + thresh_new;
% disp('halt');