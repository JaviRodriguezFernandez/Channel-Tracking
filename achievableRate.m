function [r] = achievableRate (d,K,L,Cnk,H,Ph,Wh)
r = zeros(K,1);
for k=1:K
    H_k = H{k};
    Ph_k = Ph{k};
    Wh_k = Wh{k};
    for l=1:L
        X = Cnk;
        for j=1:K
            Wh_j = Wh{j};
            if (j~=k)
                X = X + H_k(:,:,l)'*Wh_j(:,:,l)*Wh_j(:,:,l)'*H_k(:,:,l);
            end
        end
        X = Ph_k(:,:,l)'*X*Ph_k(:,:,l);
        r(k) = r(k) + log2(det(eye(d) + pinv(X)*Ph_k(:,:,l)'*H_k(:,:,l)'*Wh_k(:,:,l)*Wh_k(:,:,l)'*H_k(:,:,l)*Ph_k(:,:,l)));
    end
end
r = sum(r)/L;
end