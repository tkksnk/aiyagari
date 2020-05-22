function [r1,vmat0,mu0] = driver(r0,vmat0,mu0,Ge,Pe,knotsb,BETA,ALPHA,DELTA,SIGMA,znow,lnow,critin,critmu)

ne = size(Ge,1);
nb = size(knotsb,1);

vmat1 = zeros(nb,ne);
gmat0 = zeros(nb,ne);

% Factor Price
m0 = lnow*(r0/ALPHA/znow)^(1/(ALPHA-1));
w0 = (1-ALPHA)*znow*m0^(ALPHA)*lnow^(-ALPHA);

% instanteneous utility
umat0 = zeros(nb,nb,ne);
for ie = 1:ne
    for ib = 1:nb
        for jb = 1:nb

            enow = Ge(ie);
            know = knotsb(ib);
            kp = knotsb(jb);
            cnow = w0*enow + (1.0+r0-DELTA)*know - kp;

            if (cnow>0)
                umat0(jb,ib,ie) = log(cnow);
            else
                umat0(jb,ib,ie) = -inf;
            end

        end
    end
end

% value function iteration
diff = 1e+4;
iterin = 0;

while (diff>critin)

    for ie = 1:ne

        utemp = umat0(:,:,ie);
        vcond = Pe(ie,1)*vmat0(:,1) + Pe(ie,2)*vmat0(:,2);
        vtemp = utemp + BETA*vcond;
        [vmat1(:,ie) ivec] = max(vtemp);

%     for ie = 1:ne
%         
%         for ib = 1:nb
% 
%             utemp = umat0(:,ib,ie);
%             vcond = Pe(ie,1)*vmat0(:,1) + Pe(ie,2)*vmat0(:,2);
%             vtemp = utemp + BETA*vcond;
%             vmat1(ib,ie) = max(vtemp);
%             
%         end

    end

    diff = max(max(abs(vmat1-vmat0)));
    iterin = iterin+1;
    vmat0 = vmat1;
    % disp([iterin diff])

end

% policy function
for ie = 1:ne

    for ib = 1:nb

        utemp = umat0(:,ib,ie);
        vcond = Pe(ie,1)*vmat0(:,1) + Pe(ie,2)*vmat0(:,2);
        vtemp = utemp + BETA*vcond;
        [~,jb] = max(vtemp);
        gmat0(ib,ie) = knotsb(jb);

    end

end

% transition matrix
%AA = sparse(nb*ne,nb*ne);
AA = zeros(nb*ne,nb*ne);
wb = zeros(nb,ne);

for ie = 1:ne

    for ib = 1:nb

        know = knotsb(ib);
        kp = gmat0(ib,ie);
        kb = gridlookup2(kp,knotsb);
        wb(ib,ie) = (knotsb(kb+1)-kp)/(knotsb(kb+1)-knotsb(kb));

        for je = 1:ne

            ia = nb*(ie-1)+ib;
            ja = nb*(je-1)+kb;
            AA(ia,ja)   = wb(ib,ie)*Pe(ie,je);
            AA(ia,ja+1) = (1.0-wb(ib,ie))*Pe(ie,je);

        end

    end

end

% distribution
diffmu = 1e+4;
dist0 = reshape(mu0,2*nb,1); 

while (diffmu>critmu)

    dist1 = AA'*dist0;
    diffmu = max(abs(dist1-dist0));
    dist0 = dist1/sum(dist1);

end

mu0 = reshape(dist0,nb,2);

% Calculate K
m1 = 0.0;
for ie = 1:ne

    m1 = m1 + mu0(:,ie)'*knotsb;

end

r1 = (ALPHA)*znow*m1^(ALPHA-1)*lnow^(1-ALPHA);

end