c example for S4 iteration sweep in xy-geometry.
c Isotropic sources and specular reflection.
   do 20 i=1,iix
   do 10 j=1,jjy
10 fi(i,j) =0.0
c
   do 20 m=1,6
20 vafn(i,m) = vaf(i,m)
c
c  Start j-loop for negative eta
   do 30 jj=1,jjy   
   j = jjy1-jj
c  Start at east boundary, negative mju
   do 30 m=1,3
   a1mdy=a1*mdy(j,m)
   mdyhaf=mdy(j,m)*haf(j,m)
   call inner1(den(1,j,m),edx(1,m),vafn(1,m),q(1,j),fi(1,j))
c  Return directly from west boundary, positive mju.
   call inner2(den(1,j,m),edx(1,m),vafn(1,m+3),q(1,j),fi(1,j))
30 hafn(j,m) = mdyhaf/mdy(j,m)
c  Start j-loop for positive eta.
   do 50 j=1,jjy
c  Start at east boundary, negative mju.
   do 40 m=1,3
   a1mdy=a1*mdy(j,m)
   mdyhaf=mdy(j,m)*haf(j,m+3)
   call inner1(den(1,j,m),edx(1,m),vafn(1,m),q(1,j),fi(1,j))
C  Return directly from west boundary, positive mju.
   call inner2(den(1,j,m),edx(1,m),vafn(1,m+3),q(1,j),fi(1,j))
40 haf(j,m+3)=mdyhaf/mdy(j,m)
C Complete fluxes (all weights are equal in S4)
   do 50 i=1,iix
50 fi(i,j)=wm*fi(i,j)
C

   subroutine inner1(den, edx, vaf, q, fi)
   dimension den(1),edx(1),vaf(1),q(1),fi(1)
   common    /....../iix,iix1,.....,a1,a1a,.....,a1mdy,mdyhaf,...
C  decreasing
   i=iix1
   ii= - 1
   ilimit=1
   go to 10
!  increasing
   entry inner2(den, edx, vaf, q, fi)
   i=0
   ii=1
   ilimit=iix
10 i=i+ii
   fim=den(i)*(mdyhaf+edx(i)*vaf(i)+q(i))
   mdyhaf = a1mdy*fim-a1a*mdyhaf
   vaf(i) =  a1*fum-a1a*vaf(i)
   fi(i)= fi(i)+ fim
   if (ilimit-i) 10,20,10
20 return
   end
