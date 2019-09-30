
function Xn = ReSampleCurve(X,N)

    [n,T] = size(X);
    del(1)=0;
    for r = 2:T
        del(r) = norm(X(:,r) - X(:,r-1));
    end
    cumdel = cumsum(del)/sum(del);
    
    %keyboard;
    
    
    newdel = [1:N]/(N);
    for j=1:n
        Xn(j,:) = spline(cumdel,X(j,1:T),newdel);
    end
    q = curve_to_q(Xn);
    qn = ProjectC(q);
    Xn = q_to_curve(qn);