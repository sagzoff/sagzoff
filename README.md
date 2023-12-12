    HARMONIC ANALYSIS

x = [0 pi/3 2*pi/3 pi 4*pi/3 5*pi/3];
y = [0.8 0.6 0.4 0.7 0.9 1.1];
T = 2*pi;
w = 2*pi/T;
a0 = 2 * mean(y);
f = (0.5 * a0);
k = 2;
for n = 1:k
    an = 2 * mean(y .* cos(n * w * x)) / (0.5 * a0);
    bn = 2 * mean(y .* sin(n * w * x)) / (0.5 * a0);
    f = f + an * cos(n * w * t) + bn * sin(n * w * t);
end

vpa(f, 3);

t = linspace(x(1), x(end)); 
u = eval(f);
plot(t, u, 'r');
hold on
plot(x, y, 'go');
title('Harmonic Analysis');
legend('Harmonic Fit', 'Given Data');
xlabel('x');
ylabel('f(x)');

LANGRANGES INTERPOLATION

x=[5 6 9 11];
y=[12 13 14 16];
n=length(x);
s=0;
for i=1:n
    p=1;
    for j=1:n
        if i~=j
            p=p* (t-x(j))/(x(i)-x(j));
        end
    end
    s=s+p*y(i);
end
poly=expand(s);
fprintf('The interpolating polynomial is \t');
disp (vpa(poly, 3));
app=vpa(subs(s, t, 10),5)
t=linspace(x(1),x(end));
u=eval(poly);
plot(t,u, 'r');
hold on


NEWTON DIVIDENT 

n = length(y);
    dd(:, 1) = y;
    
    for k = 1:n
        for j = 1:n - k
            dd(j, k + 1) = (dd(j + 1, k) - dd(j, k)) / (x(j + k) - x(j));
        end
    end
    
    poly = y(1);
    
    for i = 2:n
        prod = 1;
        for j = 1:i - 1  % Fix the loop range
            prod = prod * (t - x(j));
        end
        poly = poly + prod * dd(1, i);
    end
    
    poly = simplify(poly); % Assign the result of simplify to poly
    
    t = linspace(x(1), x(end));
    z = eval(poly);
    
    
    plot(t, z);
    hold on;
    plot(x, y, 'g');  
    legend('NDD', 'Given Data');
    title('Newton''s Divided Difference Interpolation');
end

ONE THIRD RULE

% One third rule numerical integration 
f=input('Enter the function: ');
N=input('Enter the value of n :');
a=input('Enter lower limit :');
b=input('Enter upper limit :');
h=(b-a)/N;
oddsum=0;
evensum=0;
for i=1:N
    if mod(i,2)==0
        evensum=evensum+f(a+i*h);
    else
        oddsum=oddsum+f(a+i*h);
    end
end
OTRA=(h/3)*(f(a)+f(b)+2*evensum+4*oddsum); 
fprintf('\nNI by one thrid rule is %f', OTRA)

PERIODIC

%periodic functions 
x=linspace(0,2*pi); %interval
y=sin(x); % function
plot(x,y,'r')

THREE BY EIGTHT


% Newton-Cotes family of numerical integration - Code 
f=input('Enter the function: ');
N=input('Enter the value of n :');
a=input('Enter lower limit :');
b=input('Enter upper limit :');
h=(b-a)/N;
sum3=0;
sumnot3=0;
for i=1:N
    if mod(i,3)==0
        sum3=sum3+f(a+i*h);
    else
        sumnot3=sumnot3+f(a+i*h);
    end
end
TERA=(3*h/8)*(f(a)+f(b)+2*sum3+3*sumnot3); 
fprintf('\nNI by three eighth is %f', TERA)

TRAPEZODIAL

﻿% Trapezoidal Rule numerical integration 
f=input('Enter the function: ');
N=input('Enter the value of n :');
a=input('Enter lower limit :');
b=input('Enter upper limit :'); 
h=(b-a)/N;
sum=0;
for i=1:N
    sum=sum+f(a+i*h);
end
TRA=(h/2).*((2. *sum)+f(a)+f(b));
fprintf('NI by Trapezoial method is %f',TRA)

VECTOR

x=[5 6 9 11];
y=[12 13 14 16];
plot(x,y,'r')
xlabel('X')
ylabel('Y')



