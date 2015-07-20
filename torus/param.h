c      version: 9 May 2005

      parameter (jmax=512,jmax1=jmax+1,jmax2=jmax+2,jmax3=jmax-1,
     1           kmax=512,kmax1=kmax+1,kmax2=kmax+2,kmax3=kmax-1,
     2           lmax=16,lmax1=lmax/2+1,lmax2=lmax-2,
     3           max=25000,jkm1=jmax+kmax-1,neq=8)

c     32x32     wfw(452) 
c     64x64     wfw(1030)
c     128x128   wfw(2312)
c     256x256   wfw(5130)
c     512x512   wfw(11276)
c     1024x1024 wfw(24590)

