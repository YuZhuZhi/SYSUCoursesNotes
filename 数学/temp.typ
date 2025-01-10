
已知$A(2,0)$，那么设两条直线的参数方程为：

$ l_1:cases(x = 2 + l cos theta, y = l sin theta) quad, quad l_2:cases(x = 2 - l sin theta, y = l cos theta) $

代入双曲线$x^2-3y^2=3$得关于$l$的二次方程：

$ cases(l_1: l^2(c^2-3s^2) + l dot 2c + 1 = 0, l_2: l^2(s^2-3c^2) - l dot 2s + 1 = 0) $

对于二次方程$a x^2+b x+c=0$，显然有$x_1-x_2 = frac(sqrt(Delta),a)$，因此得$D E,M N$的长度：

$ cases(Delta_1=4c^2-4(c^2-3 s^2)=12s^2, Delta_2=4s^2-4(s^2-3 c^2)=12c^2) $

$ cases(|D E|=|frac(sqrt(Delta_1),c^2-3s^2)|=|frac(2 sqrt(3) s,c^2-3s^2)|, |M N|=|frac(sqrt(Delta_2),s^2-3c^2)|=|frac(2 sqrt(3) c,s^2-3c^2)|) $

所以根据不等式有$|D E|+|M N| >= 2sqrt(|D E|dot|M N|)$，计算根号内内容：

$ |D E|dot|M N| &= frac(12s c, (c^2-3s^2)(s^2-3c^2)) \ 
&= frac(12, (1/t-3t)(t-3/t)) $

计算分母内内容：

$ ... &= 1-3/t^2-3t^2+9 \ 
&= 10-(3t^2+3/t^2) \ 
&<= 10-2 sqrt(3t^2 dot 3/t^2) \ 
&= 10-6=4 $

取等条件为$3t^2=3/t^2$，即$theta=plus.minus pi/4$。分母最大值为$4$，则$|D E|dot|M N|$最小值为$3$
