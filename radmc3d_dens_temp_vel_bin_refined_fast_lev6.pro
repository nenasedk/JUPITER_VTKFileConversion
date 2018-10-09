pro radmc3d_dens_temp_vel_bin_refined_fast_lev6,nb=nb,au=au,mstar=mstar

; to get everything in CGS units from code units
CGS_temp=((1.49598d+13*au)/(sqrt((1.49598d+13*au)^3 / 6.67259d-8 / (mstar*1.9891d33))))^2/8.314d+07
CGS_dens=(mstar*1.9891d33/(1.49598d+13*au)^3)
orb_period=2*!dpi*sqrt((1.49615d+13*au)^3/(6.67259d-8*mstar*1.9891d33))
CGS_vel=1.49615d+13*au/((orb_period/(2*!dpi)))

;------------------------Numberdens:
;  divide by (mu*proton mass), then multiply the abudance of CO *1e-4 in comparison to hydrogen
; long(680.*214.*40.)+long(120.*120.*68.)+long(120.*120.*124.)+long(120.*120.*172.)+long(120.*120.*172.)+long(120.*120.*172.)+long(120.*120.*172.)

; READ in :

nbs=strcompress(string(nb),/remove_all)
print,nb

; --------------------------------- LEVEL 0
dens=dblarr(680,215,20)
temp=dblarr(680,215,20)
velo=dblarr(680,215,20,3)
dum=dblarr(long(680.*214.*40.))
dumt=dblarr(long(680.*214.*40.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_0_0.dat"
readu,1,dens
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_0_0.dat"
readu,2,temp
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_0_0.dat"
readu,3,velo
close,/all

dens=dens[*,1:214,*]
temp=temp[*,1:214,*]
velo=velo[*,1:214,*,*]
d=[[[dens]],[[reverse(dens,3)]]]*CGS_dens
t=[[[temp]],[[reverse(temp,3)]]]*CGS_temp
v=[[[velo]],[[reverse(velo,3)]]]*CGS_vel
dens=0
temp=0
velo=0

; ---------- WRITE OUT THE RADMC3D FILES --------

idum=0L
openw,101,"./output00"+nbs+"/gas_density.binp"
openw,102,"./output00"+nbs+"/dust_density.binp"
openw,103,"./output00"+nbs+"/numberdens_co.binp"
openw,104,"./output00"+nbs+"/gas_temperature.binp"
openw,105,"./output00"+nbs+"/dust_temperature.bdat"
writeu,101,1LL,8LL
writeu,101,18492800LL,1LL
writeu,102,1LL,8LL
writeu,102,18492800LL,1LL
writeu,103,1LL,8LL
writeu,103,18492800LL,1LL
writeu,104,1LL,8LL
writeu,104,18492800LL,1LL
writeu,105,1LL,8LL
writeu,105,18492800LL,1LL
for i=0,680-1 do begin
for k=0,40-1 do begin
for j=0,214-1 do begin
dum[idum]=d[i,j,k]
dumt[idum]=t[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum
writeu,102,dum/100.
writeu,103,dum*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt
writeu,105,dumt
close,/all
d=0
t=0
dum=0
dumt=0

;--------------------- VELO 0

; this is good:
;if dim eq 0 then dimr=2
;if dim eq 1 then dimr=0
;if dim eq 2 then dimr=1

dum=dblarr(long(680.*214*40.*3.))
idum=0L
openw,101,"./output00"+nbs+"/gas_velocity.binp"
writeu,101,1LL,8LL
writeu,101,18492800LL,1LL
for i=0,680-1 do begin
for k=0,40-1 do begin
for j=0,214-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum[idum]=v[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum
close,/all
v=0
dum=0

; --------------------------------- LEVEL 1
;lev 1:
dens1=dblarr(120,120,34)
temp1=dblarr(120,120,34)
velo1=dblarr(120,120,34,3)
dum1=dblarr(long(120.*120.*68.))
dumt1=dblarr(long(120.*120.*68.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_1_1.dat"
readu,1,dens1
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_1_1.dat"
readu,2,temp1
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_1_1.dat"
readu,3,velo1
close,/all

d1=[[[dens1]],[[reverse(dens1,3)]]]*CGS_dens
t1=[[[temp1]],[[reverse(temp1,3)]]]*CGS_temp
v1=[[[velo1]],[[reverse(velo1,3)]]]*CGS_vel
dens1=0
temp1=0
velo1=0

;lev1 density update
;------------------------Numberdens:
;  divide by (mu*proton mass), then multiply the abudance of CO *1e-4 in comparison to hydrogen
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,68-1 do begin
for j=0,120-1 do begin
dum1[idum]=d1[i,j,k]
dumt1[idum]=t1[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum1
writeu,102,dum1/100.
writeu,103,dum1*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt1
writeu,105,dumt1
;print,idum-1L
close,/all
d1=0
t1=0
dum1=0
dumt1=0

; ----------------velo 1 

dum1=dblarr(long(120.*120.*68.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,68-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum1[idum]=v1[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum1
close,/all
v1=0
dum1=0

; ---------------- LEVEL 2 
;lev 2:
dens2=dblarr(120,120,62)
temp2=dblarr(120,120,62)
velo2=dblarr(120,120,62,3)
dum2=dblarr(long(120.*120.*124.))
dumt2=dblarr(long(120.*120.*124.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_2_2.dat"
readu,1,dens2
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_2_2.dat"
readu,2,temp2
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_2_2.dat"
readu,3,velo2
close,/all
d2=[[[dens2]],[[reverse(dens2,3)]]]*CGS_dens
t2=[[[temp2]],[[reverse(temp2,3)]]]*CGS_temp
v2=[[[velo2]],[[reverse(velo2,3)]]]*CGS_vel
dens2=0
temp2=0
velo2=0

;lev2 density update
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,124-1 do begin
for j=0,120-1 do begin
dum2[idum]=d2[i,j,k]
dumt2[idum]=t2[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum2
writeu,102,dum2/100.
writeu,103,dum2*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt2
writeu,105,dumt2
close,/all
d2=0
t2=0
dum2=0
dumt2=0
; -----------vel 2

dum2=dblarr(long(120.*120.*124.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,124-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum2[idum]=v2[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum2
close,/all
v2=0
dum2=0

; ---------------- LEVEL 3
;lev 3:
dens3=dblarr(120,120,86)
temp3=dblarr(120,120,86)
velo3=dblarr(120,120,86,3)
dum3=dblarr(long(120.*120.*172.))
dumt3=dblarr(long(120.*120.*172.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_3_3.dat"
readu,1,dens3
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_3_3.dat"
readu,2,temp3
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_3_3.dat"
readu,3,velo3
close,/all
d3=[[[dens3]],[[reverse(dens3,3)]]]*CGS_dens
t3=[[[temp3]],[[reverse(temp3,3)]]]*CGS_temp
v3=[[[velo3]],[[reverse(velo3,3)]]]*CGS_vel
dens3=0
temp3=0
velo3=0radmc3d_dens_temp_vel_bin_refined_fast_lev6.
;lev3 density update
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
dum3[idum]=d3[i,j,k]
dumt3[idum]=t3[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum3
writeu,102,dum3/100.
writeu,103,dum3*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt3
writeu,105,dumt3
close,/all
d3=0
t3=0
dum3=0
dumt3=0
; ---------------vel 3


dum3=dblarr(long(120.*120.*172.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum3[idum]=v3[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum3
close,/all
v3=0
dum3=0

;----------------- LEVEL 4

;lev 4:
dens4=dblarr(120,120,86)
temp4=dblarr(120,120,86)
velo4=dblarr(120,120,86,3)
dum4=dblarr(long(120.*120.*172.))
dumt4=dblarr(long(120.*120.*172.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_4_4.dat"
readu,1,dens4
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_4_4.dat"
readu,2,temp4
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_4_4.dat"
readu,3,velo4
close,/all
d4=[[[dens4]],[[reverse(dens4,3)]]]*CGS_dens
t4=[[[temp4]],[[reverse(temp4,3)]]]*CGS_temp
v4=[[[velo4]],[[reverse(velo4,3)]]]*CGS_vel
dens4=0
temp4=0
velo4=0
;lev4 density update
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
dum4[idum]=d4[i,j,k]
dumt4[idum]=t4[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum4
writeu,102,dum4/100.
writeu,103,dum4*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt4
writeu,105,dumt4
close,/all
d4=0
t4=0
dum4=0
dumt4=0

;----------------------------- VEL 4
dum4=dblarr(long(120.*120.*172.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum4[idum]=v4[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum4
close,/all
v4=0
dum4=0

; ----------------------------- LEVEL 5 

;lev 5:
dens5=dblarr(120,120,86)
temp5=dblarr(120,120,86)
velo5=dblarr(120,120,86,3)
dum5=dblarr(long(120.*120.*172.))
dumt5=dblarr(long(120.*120.*172.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_5_5.dat"
readu,1,dens5
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_5_5.dat"
readu,2,temp5
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_5_5.dat"
readu,3,velo5
close,/all
d5=[[[dens5]],[[reverse(dens5,3)]]]*CGS_dens
t5=[[[temp5]],[[reverse(temp5,3)]]]*CGS_temp
v5=[[[velo5]],[[reverse(velo5,3)]]]*CGS_vel
dens5=0
temp5=0
velo5=0

;lev5 density update
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
dum5[idum]=d5[i,j,k]
dumt5[idum]=t5[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum5
writeu,102,dum5/100.
writeu,103,dum5*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt5
writeu,105,dumt5
close,/all
d5=0
t5=0
dum5=0
dumt5=0

;----------------------------- VEL 5
dum5=dblarr(long(120.*120.*172.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum5[idum]=v5[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum5
close,/all
v5=0
dum5=0

; ---------------- LEVEL 6 
;lev 6:
dens6=dblarr(120,120,86)
temp6=dblarr(120,120,86)
velo6=dblarr(120,120,86,3)
dum6=dblarr(long(120.*120.*172.))
dumt6=dblarr(long(120.*120.*172.))
openr,1,"./output00"+nbs+"/gasdensity"+nbs+"_6_6.dat"
readu,1,dens6
openr,2,"./output00"+nbs+"/gastemperature"+nbs+"_6_6.dat"
readu,2,temp6
openr,3,"./output00"+nbs+"/gasvelocity"+nbs+"_6_6.dat"
readu,3,velo6
close,/all
d6=[[[dens6]],[[reverse(dens6,3)]]]*CGS_dens
t6=[[[temp6]],[[reverse(temp6,3)]]]*CGS_temp
v6=[[[velo6]],[[reverse(velo6,3)]]]*CGS_vel
dens6=0
temp6=0
velo6=0

;lev6 density update
idum=0L
openu,101,"./output00"+nbs+"/gas_density.binp",/append
openu,102,"./output00"+nbs+"/dust_density.binp",/append
openu,103,"./output00"+nbs+"/numberdens_co.binp",/append
openu,104,"./output00"+nbs+"/gas_temperature.binp",/append
openu,105,"./output00"+nbs+"/dust_temperature.bdat",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
dum6[idum]=d6[i,j,k]
dumt6[idum]=t6[i,j,k]
idum=idum+1L
endfor
endfor
endfor
writeu,101,dum6
writeu,102,dum6/100.
writeu,103,dum6*1.e-4/(2.3*1.6726e-24)
writeu,104,dumt6
writeu,105,dumt6
close,/all
d6=0
t6=0
dum6=0
dumt6=0
;----------------------------- VEL 6

dum6=dblarr(long(120.*120.*172.*3.))
idum=0L
openu,101,"./output00"+nbs+"/gas_velocity.binp",/append
for i=0,120-1 do begin
for k=0,172-1 do begin
for j=0,120-1 do begin
for dim=0,2 do begin
if dim eq 0 then dimr=2
if dim eq 1 then dimr=0
if dim eq 2 then dimr=1
dum6[idum]=v6[i,j,k,dimr]
idum=idum+1L
endfor
endfor
endfor
endfor
writeu,101,dum6
close,/all
v6=0
dum6=0


end
