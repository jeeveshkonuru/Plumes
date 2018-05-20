function [z,b,u,rho,x] = plumes_makeprof(sout)
%THIS ARGUMENTS ARE USED IN ACTUAL FUNCTION
if sout.nondim == 0
    [z,b,u,rho,x] = plumes_run(sout);
elseif sout.nondim == 1
    [z,b,u,rho,x] = plumes_run_nondim(sout);
end