//Date : 27/01/2022
//Aim : To solve a second order ordinary differential equation with general boundary conditions using finite difference method
clc
clear
clf
//y''-y=0 with y(0)=0 and y(2)=3.62686        DD
//y''-y=0 with y(0)=1 and y'(1)=0             DN
//y''+y=sin(3x) with y'(0)=0 and y(pi)=0.5      ND
//y''+y=sin(3x) with y'(0)=0 and y'(pi)=1      NN
//y''+y=sin(3x) with y(0)+y'(0)=-1 and y'(pi)=1   Robin
function cy_prime = f(x)
    cy_prime=0
endfunction
function cy = g(x)
    cy=1
endfunction
function R = r(x)
    R=sin(3*x)
endfunction
funcprot(0)
x0=input("Initial value of x : ")
xn=input("Boundary value of x : ")
n=10
//h=0.1//step size
h=(xn-x0)/n
x=x0:h:xn
//n=int((xn-x0)/h)
disp("n : ",n)//number of intervals 
//define general boundary conditions
alpha1=input("The value of alpha 1 :")
alpha2=input("The value of alpha 2 :")
alpha3=input("The value of alpha 3 :")
Beta1=input("The value of Beta 1 :")
Beta2=input("The value of Beta 2 :")
Beta3=input("The value of Beta 3 :")
A=zeros(n+1,n+1)//coefficient matrix
B=zeros(n+1,1)//column matrix
for i=1:n+1
    if alpha1==0 && Beta1==0 then      //NN
        if i==1
            A(i,i)=(-2/(h^2))+g(x(i))
            A(i,i+1)=2/(h^2)
            B(i)=r(x(i))+(((2/h)-f(x(i)))*(alpha3/alpha2))
        elseif i==n+1
            A(i,i)=(-2/(h^2))+g(x(i))
            A(i,i-1)=2/(h^2)
            B(i)=r(x(i))-(((2/h)+f(x(i)))*(Beta3/Beta2))
        else
            A(i,i)=(-2/h^2)+g(x(i))
            A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
            A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
            B(i)=r(x(i))
        end
    elseif alpha2==0 && Beta2==0        //DD
        if i==1 then
            A(i,i)=1
            B(i)=alpha3
        elseif i==n+1
            A(i,i)=1
            B(i)=Beta3
        else
            A(i,i)=(-2/h^2)+g(x(i))
            A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
            A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
            B(i)=r(x(i))
        end
    elseif alpha1==0 && Beta2==0      //ND
        if i==1 then
            A(i,i)=(-2/(h^2))+g(x(i))
            A(i,i+1)=2/(h^2)
            B(i)=r(x(i))+(((2/h)-f(x(i)))*(alpha3/alpha2))
        elseif i==n+1
            A(i,i)=1
            B(i)=Beta3
        else
            A(i,i)=(-2/h^2)+g(x(i))
            A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
            A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
            B(i)=r(x(i))
        end
    elseif alpha2==0 && Beta1==0       //DN
        if i==1 then
           A(i,i)=1
           B(i)=alpha3 
        elseif i==n+1
            A(i,i)=(-2/(h^2))+g(x(i))
            A(i,i-1)=2/(h^2)
            B(i)=r(x(i))-(((2/h)+f(x(i)))*(Beta3/Beta2))
        else
            A(i,i)=(-2/h^2)+g(x(i))
            A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
            A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
            B(i)=r(x(i))
        end
    else                       //General Robin
        if i==1 then
            A(i,i)=(-2/(h^2))+g(x(i))+(((2/h)-f(x(i)))*(alpha1/alpha2))
            A(i,i+1)=2/(h^2)
            B(i)=r(x(i))+(((2/h)-f(x(i)))*(alpha3/alpha2))
        elseif i==n+1
            A(i,i)=(-2/(h^2))+g(x(i))+(((2/h)+f(x(i)))*(Beta1/Beta2))
            A(i,i-1)=2/(h^2)
            B(i)=r(x(i))-(((2/h)+f(x(i)))*(Beta3/Beta2))
       else
           A(i,i)=(-2/h^2)+g(x(i))
           A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
           A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
           B(i)=r(x(i))
        end
    end
end
disp("The coefficient matrix on LHS : ",A)
disp("The column matrix on RHS : ",B)
y=inv(A)*B
disp("The solution of the given BVP is : ",y)
plot(x,y')
//exact solution
for j=1:n+1       //for ND
    solution(j) = (-cos(x(j)))+((3/8)*sin(x(j)))-((1/8)*sin(3*x(j)))
end
disp("Exact solution : ",solution)
plot(x,solution','*r')
title("Solving Boundary valued Problem using Finite Difference Method")
title color Red
title fontsize 4
xlabel("X----->")
xlabel color magenta fontsize 4
ylabel("Y(X)---->")
ylabel color magenta fontsize 4
legend(["Using Finite Difference Method";"Exact Solution"])
xgrid
