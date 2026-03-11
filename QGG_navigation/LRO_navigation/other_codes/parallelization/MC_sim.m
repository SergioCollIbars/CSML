clear; clc;
MC = 100;
pool = gcp;

metaKernelPath = '/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm';

addpath("data/")
addpath(genpath("functions/"))

option = "SQ"; % options: parallel (PL) or sequential (SQ)

if option == "PL"
    tic
    futures(MC) = parallel.FevalFuture;
    for k = 1:MC
        futures(k) = parfeval(pool, @run_one_MC_spice, 1, metaKernelPath);
    end
    
    results = cell(MC,1);
    for n = 1:MC
        [idx, out] = fetchNext(futures);
        results{idx} = out;
        fprintf("Done %d/%d (run %d)\n", n, MC, idx);
    end
    toc
elseif(option == "SQ")
    tic
    for k =1:MC
        run_one_MC_spice(metaKernelPath);
    end
    toc
end