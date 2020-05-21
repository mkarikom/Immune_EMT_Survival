function [scaled,scaleoffset] = findoffset(x)
    n = length(x(:));
    allequal = true;
    offprobs = zeros(size(x));
    nsig = 1;
    while allequal
        nsig = nsig * 10;
        nextsig = nsig * 10;
        scaleoffset = ceil(x(1,1)*nsig)/nsig;
        offset = ones(size(x))*scaleoffset;
        nextoffset = ones(size(x))*ceil(x(1,1)*nextsig)/nextsig;
        offprobs = nextoffset-offset;
        sigoffprobs = ceil(offprobs*nextsig)/nextsig;
        if length(unique(sigoffprobs)) > 1
            allequal = false;
        end
        scaled=x-offset;
    end
end