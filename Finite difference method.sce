//Date :20/01/2022
//Aim : To solve a second order ordinary differential equation with Dirichlet boundary conditions using Finite Difference method
clc
clear
function Cy_prime = f(x)
    Cy_prime = 0
endfunction
function Cy = g(x)
    Cy = -1
endfunction
function R = r(x)
    R = 0
endfunction
funcprot(0)
x0=input("First initial value of x : ")
xn=input("Boundary value of x : ")
h=0.05
x=x0:h:xn
n=(xn-x0)/h
//disp(n)
A =zeros (n+1,n+1)//coefficient matrix
for i=1:n+1//calculating the coefficient matrix
    if i==1 || i==n+1
        A(i,i)=1
    else
        A(i,i)=(-2/h^2)+g(x(i))
        A(i,i-1)=(1/h^2)-(f(x(i))/(2*h))
        A(i,i+1)=(1/h^2)+(f(x(i))/(2*h))
end
end
disp("Coefficient matrix : ",A)
B=zeros(n+1,1)
y0=input("First initial value of y : ")
yn=input("Boundary value of y : ")
for j = 1:n+1//calculating B matrix
    if j==1
        B(j,1)=y0
    elseif j==n+1
        B(j,1)=yn
    else
        B(j,1)=r(x(j))
end
end
disp("B matrix : ",B)
solution = inv(A)*B
disp("The solution : ",solution)
plot(x',solution,'*r')
for k=1:n+1
    exact(k)=(0.5*exp(x(k)))-(0.5*exp(-x(k)))
    disp(exact(k))
end
plot(x',exact)
title("Solving Boundary valued Problem using Finite Difference Method")
title color Red
title fontsize 4
xlabel("X----->")
xlabel color magenta fontsize 4
ylabel("Y(X)---->")
ylabel color magenta fontsize 4
legend(["Using Finite Difference Method";"Exact Solution"])
xgrid

