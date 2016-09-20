f1(x)=a*1/x+b*1/x**2+c*1/x**3 +d*1/x**4 +e*1/x**5+f*x+g
	a=1; b=2; c=3; d=2;e=1; f=1;g=1
#f2(x)=a2*x+b2*x**2+c2*x**3 +d2*x**4 +e2*x**5+f2
#	a2=1; b2=2; c2=3; d2=2; e2=1; g2=1; f2=1; 
N=35
neff=N-(3.1311804+0.1784/(N-3.1311804)**2)
fit f1(x) 'phase_shifts_RB.txt' using (2/($1**2+1/neff**2)):(($1**3)/(-tan($6))) via a,b ,c ,d ,e,f,g
#fit f1(x) 'phase_shifts_RB.txt' using 1:(-($1**3)/tan($5)) via a,b,c  ,d #,e,f
#fit f2(x) 'phase_shifts_RB.txt' using 1:(-($1**3)/tan($6)) via a2,b2,c2 ,d2 ,e2, f2

set terminal epslatex color font 'Helvetica,10'
#set output "cali.eps"
set title "Calibration Curve" font ",25"     # note newline!
set lmargin at screen 0.08
set rmargin at screen 0.95
set bmargin at screen 0.10
set tmargin at screen 0.90
set pointsize 1.5                              # larger point!
set xlabel  '$R[a_0]$' font ",20" offset 0,0.5          # Greek symbols!
set ylabel  'Scattering Length $[a_0]$'  font ",20" offset 0.5,0         # italics!
set tics font ", 15"
set key right spacing 2.0 font ",18"
plot 1/f1(x) title 'fitting of P-Triplet',\
     'phase_shifts_RB.txt' using (2/($1**2+1/neff**2)):(tan($6)/(-($1**3))) title 'P-Triplet'
     #'phase_shifts_RB.txt' using 1:(-tan($5)/($1**3)) title 'P-Singulet'
#     'phase_shifts_RB.txt' using 1:(-tan($6)/($1**3)) title 'P-Triplet'
#f1(x) title 'fitting of P-Triplet',\
#     'phase_shifts_RB.txt' using 1:(-tan($5)/($1**3)) title 'P-Triplet'
