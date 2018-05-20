function S=param_file_to_struct(paramfile)
    eval(paramfile);
    S=sout;
end