set term x11
set macros

# name of files in ./hubTri2Dcondres/ to plot
filen="n08u120tp00_ndk002"


n0=system("perl -e \'$a=\"/".filen."\"; $a=~s/^.n//; $a=~s/u.+$//; print \"$a\"; exit\'")
n0=n0*1
print "Number of sites n0= ",n0


##################
print "ploting nf"

set pm3d map
set cont base
set cntrparam levels incremental 0, 1, 40
unset surface
set table 'cont.dat'
splot "./hubTri2Dcondres/".filen."_nf.dat" u 1:2:3 w l
unset table
#pause -1
reset

set xlabel "mu [t]"
set ylabel "T [t]"
set title "fermionic density/number plot "
set pm3d map
splot "./hubTri2Dcondres/".filen."_nf.dat" u 1:2:($3)  with pm3d ,\
'cont.dat' w l t "contours",\
"./hubTri2Dcondres/".filen."_nfd100" u 3:2:(0.0) w l t "n=1. (half-filling)",\
"./hubTri2Dcondres/".filen."_nfd101" u 3:2:(0.0) w l t "n=0.95",\
"./hubTri2Dcondres/".filen."_nfd102" u 3:2:(0.0) w l t "n=0.90",\
"./hubTri2Dcondres/".filen."_nfd103" u 3:2:(0.0) w l t "n=0.85",\
"./hubTri2Dcondres/".filen."_nfd104" u 3:2:(0.0) w l t "n=0.80",\
"./hubTri2Dcondres/".filen."_nfd105" u 3:2:(0.0) w l t "n=0.75",\
"./hubTri2Dcondres/".filen."_nfd106" u 3:2:(0.0) w l t "n=0.70" 
pause -1



##################
print "Ploting tau"
set xlabel "T [t]"
set ylabel "tau, sum rule [t]"
set yrange[0:]

plot \
 "./hubTri2Dcondres/".filen."_tau_nfd100" u 2:($4)  w lp t "n=1.00",\
 "./hubTri2Dcondres/".filen."_tau_nfd101" u 2:($4)  w lp t "n=0.95",\
 "./hubTri2Dcondres/".filen."_tau_nfd102" u 2:($4)  w lp t "n=0.90",\
 "./hubTri2Dcondres/".filen."_tau_nfd103" u 2:($4)  w lp t "n=0.85",\
 "./hubTri2Dcondres/".filen."_tau_nfd104" u 2:($4)  w lp t "n=0.80",\
 "./hubTri2Dcondres/".filen."_tau_nfd105" u 2:($4)  w lp t "n=0.75",\
 "./hubTri2Dcondres/".filen."_tau_nfd106" u 2:($4)  w lp t "n=0.70"
pause -1



##################
print "ploting optical conductivity"
#set which fermion density to plot
nnfd="nfd100"
nnfd="nfd101"
nnfd="nfd102"
nnfd="nfd103"
nnfd="nfd104"
nnfd="nfd105"
nnfd="nfd106"


nnfd="nfd103"
print "nfd = ",nnfd
set title filen." ".nnfd

set autoscale xy
set xlabel "omega [t]"
set ylabel "optical conductivity"
v1=NaN
v2=NaN
fsetv1(x)=(v1=x,v1)
fsetv2(x)=(v2=x,v2)
#for number of indices in file : STATS_blocks
stats "./hubTri2Dcondres/".filen."_condeta_".nnfd
do for[i=0:STATS_blocks-2] {
plot \
"./hubTri2Dcondres/".filen."_condeta_".nnfd index i u (fsetv2($2)*NaN):(fsetv1($1)*NaN)  w l not "index".gprintf("%g",i) ,\
"./hubTri2Dcondres/".filen."_condeta_".nnfd index i u 3:($4/n0)  w l t "T index ".gprintf("%g",i).", T=".gprintf("%g",v1).", pure j-j corelation" ,\
"./hubTri2Dcondres/".filen."_condeta_".nnfd index i u 3:($5/n0)  w l t "T index ".gprintf("%g",i).", T=".gprintf("%g",v1).", D_c(om=0) from tau"
pause -1
}





#set yrange[-1:5]
#set autoscale y
print "Ploting tau"
set xlabel "T [t]"
set ylabel "tau, sum rule [t]"


plot \
 "./hubTri2Dcondres/".filen."_tau_nfd100" u 2:($4)  w lp t "tau",\
 "./hubTri2Dcondres/".filen."_tau_nfd100" u ($2):($5)  w lp t "j-j corr spectral sum" 
pause -1


plot \
 "./hubTri2Dcondres/".filen."_tau_nfd100" u ($2):(2*$13/$4) w lp t "p=0.00 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd101" u ($2):(2*$13/$4) w lp t "p=0.05 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd102" u ($2):(2*$13/$4) w lp t "p=0.10 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd103" u ($2):(2*$13/$4) w lp t "p=0.15 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd104" u ($2):(2*$13/$4) w lp t "p=0.20 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd105" u ($2):(2*$13/$4) w lp t "p=0.25 2*D_c phase/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd106" u ($2):(2*$13/$4) w lp t "p=0.30 2*D_c phase/tau" 
pause -1

plot \
 "./hubTri2Dcondres/".filen."_tau_nfd100" u ($2):(2*$8/$4) w lp t "p=0.00 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd101" u ($2):(2*$8/$4) w lp t "p=0.05 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd102" u ($2):(2*$8/$4) w lp t "p=0.10 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd103" u ($2):(2*$8/$4) w lp t "p=0.15 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd104" u ($2):(2*$8/$4) w lp t "p=0.20 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd105" u ($2):(2*$8/$4) w lp t "p=0.25 2*D_c tau/tau" ,\
 "./hubTri2Dcondres/".filen."_tau_nfd106" u ($2):(2*$8/$4) w lp t "p=0.30 2*D_c tau/tau" 
pause -1


reset
print "ploting conduc. channels"
set xlabel "omega [t]"
set ylabel "channel weight"

fsetv1(x)=(v1=x,v1)
fsetv2(x)=(v2=x,v2)
#for number of indices in file : STATS_blocks
stats "./hubTri2Dcondres/".filen."_condchan_".nnfd 
do for[i=0:STATS_blocks-2] {
plot \
"./hubTri2Dcondres/".filen."_condchan_".nnfd index i u (fsetv2($2)*NaN):(fsetv1($1)*NaN)  w l not "index".gprintf("%g",i) ,\
"./hubTri2Dcondres/".filen."_condchan_".nnfd index i u 3:4  w l t "T index".gprintf("%g",i).", T=".gprintf("%g",v1).", pure j-j corr." ,\
"./hubTri2Dcondres/".filen."_condchan_".nnfd index i u 3:5  w l t "T index".gprintf("%g",i).", T=".gprintf("%g",v1).", D_c from tau at om=0" 
pause -1
}


exit


