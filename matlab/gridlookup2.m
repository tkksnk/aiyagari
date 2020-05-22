function ix = gridlookup2(x0,xgrid)

nx = size(xgrid,1);

ix = 0;
for jx = 1:nx
    if (x0<=xgrid(jx)); break; end;
    ix = ix+1;
end

ix = min(max(1,ix),nx-1);