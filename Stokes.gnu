t   = time(0.0)
ξ1  = 0.9                  # laser linear polarization horizontal (ξ1=1) or vertical ξ1=-1
ξ2  = 0.0                  # laser linear polarization φ=π/4 (ξ2=1) or φ=-π/4 (ξ2=-1)
ξ3  = sqrt(1-ξ1*ξ1-ξ2*ξ2)  # laser circular polarization degree ξ3 [-1:1]

Ex  = sqrt((1+ξ1)/2)
Ey  = sqrt((1-ξ1)/2)
δ   = Ex*Ey!=0 ? atan2(ξ3,ξ2) : 0
x   = Ex*sin(t)
y   = Ey*sin(t+δ)

set xrange[-1.1:1.1]; set xlabel 'x'
set yrange[-1.1:1.1]; set ylabel 'y'
set arrow 1 from 0,0 to x,y filled lw 2 
set arrow 2 from 0,0 to x,0 filled 
set arrow 3 from 0,0 to 0,y filled 
set arrow 4 from x,0 to x,y nohead dt 2
set arrow 5 from 0,y to x,y nohead dt 2
set grid
plot NaN notitle
pause 0.05
reread

