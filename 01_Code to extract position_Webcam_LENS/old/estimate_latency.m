function latencies = estimate_latency(filename)

data=importdata(filename);
v_rest_mean = mean(data(rest_indexes,9));
v_rest_std = std(data(rest_indexes,9));

thres = v_rest_mean-4*v_rest_std;
peaks = find(data(:,9)<thres);
