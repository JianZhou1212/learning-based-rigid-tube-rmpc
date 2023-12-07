function Fs = compute_mrpi_set(Ak, W, epsilon)

[nx, ~] = size(Ak); 
s = 0; 
alpha = 1000;
Ms = 1000;
it = 0;
while(alpha > epsilon/(epsilon + Ms))
    s = s + 1;
    alpha = max(W.support(Ak^s*(W.A)')./W.b);
    mss = zeros(2*nx, 1);
    for i = 1:s
        mss = mss + W.support([Ak^i, - Ak^i]);
    end
    Ms = max(mss);
    it = it+1;
end

Sam_num=10000;
Sample=50*rand(2,Sam_num)-25;


yalmip('clear')
z = sdpvar(nx,s); 
proj_sample= sdpvar(nx,1);
Input_sample=sdpvar(nx,1);
cns=[];
item=zeros(nx,1);
for i=1:s
    item=item+Ak^(i-1)*z(:,i);
    cns=[cns,W.A*z(:,i)<=W.b];
end
cns=[cns,proj_sample==item/(1-alpha)];
obj=norm(proj_sample-Input_sample);
ops = sdpsettings('relax',0);
MRPI = optimizer(cns,obj,ops,Input_sample,proj_sample);

parfor kkk=1:Sam_num
     [sample_proj(:,kkk),errorcode]=MRPI(Sample(:,kkk));
end

[convhull_index,~] = convhull(sample_proj');
MRPI_W=sample_proj(:,convhull_index);
Fs=Polyhedron(MRPI_W');

end
