function f = w_generator(varargin)
%w_generator(method,mu,Sigma) (hack) for unifrnd
    if nargin == 0
        mu = [0 0];
        Sigma = 0.1*[0.1 0; 0 1];
        f = mvnrnd(mu,Sigma,1);
        %%%DEBUG%%%
        %f = unifrnd(-0.2,0.2,[1 2]);
    else
        method = varargin{1};
        if method == "mvnrnd"
            mu = varargin{2};
            Sigma = varargin{3};
            f=mvnrnd(mu,Sigma,1);
        elseif method == "uniform"
            lb = varargin{2};
            ub = varargin{3};
            shape = varargin{4};
            f = unifrnd(lb,ub,shape);
        elseif method == "rmvnrnd"
            %rmvnrnd(MU,SIG,N,A,B)
            mu = varargin{2};
            Sigma = varargin{3};
            N = varargin{4};
            A = varargin{5};
            b = varargin{6};  
            f = rmvnrnd(mu,Sigma,N,A,b);
        else
            display("Method name isn't valid");
            abc = input("Press Enter to continue anyway.");
            f = -1;
        end
    end
end