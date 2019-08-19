!BL10
    subroutine Sweep2D(BC1,BC2,BC3,BC4,totNFM_x,totnfm_y,ngauss,ng,&
                       Nmat,SigT,Q_nij,mu,eta,w,dx,dy,fmmid2D,phi_ij)
    implicit none
    ! Variable Globale
    CHARACTER(50), intent(in) :: BC1,BC2,BC3,BC4
    integer(kind=4), intent(in) :: totNFM_x,totnfm_y,ngauss,ng,Nmat
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2,ng), intent(in) :: Q_nij 
    real(kind=8), dimension(ngauss*(ngauss+2)/8), intent(in) :: mu, eta, w
    real(kind=8), dimension(totNFM_x), intent(in) :: dx
    real(kind=8), dimension(totNFM_y), intent(in) :: dy
    real(kind=8), dimension(Nmat,ng),  intent(in) :: SigT
    integer(kind=4), dimension(totNFM_x,totnfm_y), intent(in) :: fmmid2D
    real(kind=8), dimension(totNFM_x,totnfm_y,ng), intent(out) :: phi_ij
    ! Variable locale
    real(kind=8), dimension(totNFM_x,totnfm_y,ngauss*(ngauss+2)/2,ng) :: phi_ijn,phi_aijn,phi_bijn
    real(kind=8), dimension(totNFM_x,totnfm_y,ng) :: phi
    real(kind=8), dimension(totNFM_x+1,totnfm_y,ngauss*(ngauss+2)/2) :: phix,phiax,phibx
    real(kind=8), dimension(totNFM_x,totnfm_y+1,ngauss*(ngauss+2)/2) :: phiy,phiay,phiby
    real(kind=8) :: ax,ay,nume,deno,beta
    integer(kind=4) :: i,j,n,k,numord,nn
    numord = ngauss*(ngauss+2)/2
    phi = 0.0d0
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    nn=1
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          phiax(1,j,n)   = 0.0D0 ! pour mu>0
          phibx(1,j,n)   = 1.0D0
          do i = 1,totNFM_x         
             phiay(i,1,n)   = 0.0D0
             phiby(i,1,n)   = 1.0D0
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_nij(i,j,n,k) + ax*phiax(i,j,n) + ay*phiay(i,j,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0d0
             !> end
             phiax(i+1,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i,j,n)
             !> simple negative flux fixup         
             if (phiax(i+1,j,n)<0) phiax(i+1,j,n) = 0.0d0
             !> end
             nume = Q_nij(i,j,n,k) + ax*phibx(i,j,n) + ay*phiby(i,j,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0d0
             !> end
             phibx(i+1,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i,j,n)
             !> simple negative flux fixup         
             if (phibx(i+1,j,n)<0) phibx(i+1,j,n) = 0.0d0
             !> end
             !> Boundary condition in the left side
             if (BC3 == 'Vacuum Left') phix(i,1,n) = 0.0D0
             if (BC3 == 'Reflective Left') phix(i,1,n)=phiax(i,1,n)/(1.0D0+(phiax(i,1,n)-phibx(i,1,n)))
             !> end
          enddo
          do i = 1,totNFM_x
             phiay(i,j+1,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j,n)
             phiby(i,j+1,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j,n)
             !> simple negative flux fixup         
             if (phiay(i,j+1,n)<0) phiay(i,j+1,n) = 0.0d0
             if (phiby(i,j+1,n)<0) phiby(i,j+1,n) = 0.0d0
             !> end
             !> Boundary condition in the bottom side
             if (BC2 == 'Vacuum Bottom') phiy(1,j,n) = 0.0D0
             if (BC2 == 'Reflective Bottom') phiy(1,j,n)=phiay(1,j,n)/(1.0D0+(phiay(1,j,n)-phiby(1,j,n)))
             !> end
          enddo
       enddo
    nn=nn+1
    enddo
    nn=1
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          do i = 1,totNFM_x
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             !> end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             !> simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             !> simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo

    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do j = totNFM_y,1,-1
          phiax(totNFM_x+1,j,n)  = 0.0D0 
          phibx(totNFM_x+1,j,n)  = 1.0D0
          do i = totNFM_x,1,-1         
             phiay(i,totNFM_y+1,n)  = 0.0D0
             phiby(i,totNFM_y+1,n)  = 1.0D0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_nij(i,j,n,k) + ax*phiax(i+1,j,n) + ay*phiay(i,j+1,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0d0
             !> end
             phiax(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i+1,j,n)
             !> simple negative flux fixup         
             if (phiax(i,j,n)<0) phiax(i,j,n) = 0.0d0
             !> end
             nume = Q_nij(i,j,n,k) + ax*phibx(i+1,j,n) + ay*phiby(i,j+1,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0d0
             !> end
             phibx(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i+1,j,n)
             !> simple negative flux fixup         
             if (phibx(i,j,n)<0) phibx(i,j,n) = 0.0d0
             !> end
             !> Boundary condition in the right side
             if (BC4 == 'Vacuum Right') phix(totNFM_x+1,1,n)  = 0.0D0
             if (BC4 == 'Reflective Right') phix(i+1,j,n) =phiax(totNFM_x+1,1,n)/&
                (1.0D0+(phiax(totNFM_x+1,1,n)-phibx(totNFM_x+1,1,n)))
             !> end
          enddo
          do i = totNFM_x,1,-1         
             phiay(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j+1,n)
             phiby(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j+1,n)
             !> simple negative flux fixup         
             if (phiay(i,j,n)<0) phiay(i,j,n) = 0.0d0
             if (phiby(i,j,n)<0) phiby(i,j,n) = 0.0d0
             !> end
             !> Boundary condition in the top side
             if (BC1 == 'Vacuum Top') phiy(1,totNFM_y+1,n) = 0.0D0
             if (BC1 == 'Reflective Top') phiy(i,totNFM_y+1,n)=phiay(i,totNFM_y+1,n)/&
                 (1.0D0+(phiay(1,totNFM_y+1,n)-phiby(1,totNFM_y+1,n)))
             !> end
          enddo
       enddo
    enddo
    
    do n = 1,numord/4
       do j = totNFM_y,1,-1         
          do i = totNFM_x,1,-1         
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             !> end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             !> simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1         
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             !> simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             !> end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    nn=1
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1
          phiax(1,j,n) = 0.0D0 ! pour mu>0
          phibx(1,j,n)   = 1.0D0
          do i = 1,totNFM_x         
             phiay(i,totNFM_y+1,n)   = 0.0D0
             phiby(i,totNFM_y+1,n)   = 1.0D0
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_nij(i,j,n,k) + ax*phiax(i,j,n) + ay*phiay(i,j+1,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0d0
             !> end
             phiax(i+1,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i,j,n)
             !> simple negative flux fixup         
             if (phiax(i+1,j,n)<0) phiax(i+1,j,n) = 0.0d0
             !> end
             nume = Q_nij(i,j,n,k) + ax*phibx(i,j,n) + ay*phiby(i,j+1,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0d0
             !> end
             phibx(i+1,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i,j,n)
             !> simple negative flux fixup         
             if (phibx(i+1,j,n)<0) phibx(i+1,j,n) = 0.0d0
             ! end
             !> Boundary condition in the left side
             if (BC3 == 'Vacuum Left') phix(i,j,n) = 0.0D0
             if (BC3 == 'Reflective Left') beta = 1.0D0
             !> end
             phix(i,j,n) = beta*phiax(i,j,n)/(1.0D0 + beta*(phiax(i,j,n)-phibx(i,j,n)))
          enddo
          do i = 1,totNFM_x
             phiay(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j+1,n)
             phiby(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j+1,n)
             !> simple negative flux fixup         
             if (phiay(i,j,n)<0) phiay(i,j,n) = 0.0d0
             if (phiby(i,j,n)<0) phiby(i,j,n) = 0.0d0
             !> end
             !> Boundary condition in the left side
             if (BC1 == 'Vacuum Top') beta = 0.0D0
             if (BC1 == 'Reflective Top') beta = 1.0D0
             !> end
             phiy(i,j+1,n) = beta*phiay(i,j+1,n)/(1.0D0 + beta*(phiay(i,j+1,n)-phiby(i,j+1,n)))
          enddo
       enddo
       nn=nn+1
    enddo

    nn=1
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1         
          do i = 1,totNFM_x
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             !> end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             !> simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             !> simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    nn=1
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y
          phiax(totNFM_x+1,j,n)   = 0.0D0 
          phibx(totNFM_x+1,j,n)   = 1.0D0
          do i = totNFM_x,1,-1  
             phiby(i,1,n)   = 1.0D0     
             phiay(i,1,n)   = 0.0D0
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay
             nume = Q_nij(i,j,n,k) + ax*phiax(i+1,j,n) + ay*phiay(i,j,n)
             phi_aijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_aijn(i,j,n,k)<0) phi_aijn(i,j,n,k) = 0.0d0
             !> end
             phiax(i,j,n) = 2.0D0*phi_aijn(i,j,n,k)-phiax(i+1,j,n)
             !> simple negative flux fixup         
             if (phiax(i,j,n)<0) phiax(i+1,j,n) = 0.0d0
             !> end
             nume = Q_nij(i,j,n,k) + ax*phibx(i+1,j,n) + ay*phiby(i,j,n)
             phi_bijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_bijn(i,j,n,k)<0) phi_bijn(i,j,n,k) = 0.0d0
             !> end
             phibx(i,j,n) = 2.0D0*phi_bijn(i,j,n,k)-phibx(i+1,j,n)
             !> simple negative flux fixup         
             if (phibx(i,j,n)<0) phibx(i,j,n) = 0.0d0
             !> end
             !> Boundary condition in the right side
             if (BC4 == 'Vacuum Right') beta = 0.0D0
             if (BC4 == 'Reflective Right') beta = 1.0D0
             !> end
             phix(i+1,j,n) = beta*phiax(i+1,j,n)/(1.0D0 + beta*(phiax(i+1,j,n)-phibx(i+1,j,n)))
          enddo
          do i = totNFM_x,1,-1         
             phiay(i,j+1,n) = 2.0D0*phi_aijn(i,j,n,k)-phiay(i,j,n)
             phiby(i,j+1,n) = 2.0D0*phi_bijn(i,j,n,k)-phiby(i,j,n)
             !> simple negative flux fixup         
             if (phiay(i,j+1,n)<0) phiay(i,j+1,n) = 0.0d0
             if (phiby(i,j+1,n)<0) phiby(i,j+1,n) = 0.0d0
             !> end
             !> Boundary condition in the bottom side
             if (BC2 == 'Vacuum Bottom') beta = 0.0D0
             if (BC2 == 'Reflective Bottom') beta = 1.0D0
             !>
             phiy(i,j,n) = beta*phiay(i,j,n)/(1.0D0 + beta*(phiay(i,j,n)-phiby(i,j,n)))
          enddo
       enddo
       nn=nn+1
    enddo

    nn=1
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y        
          do i = totNFM_x,1,-1         
             ax = 2.0d0*abs(mu(nn))/dx(i)
             ay = 2.0d0*abs(eta(nn))/dy(j)
             deno = SigT(fmmid2D(i,j),k) + ax + ay 
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             phi_ijn(i,j,n,k) = nume/deno
             !> simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             !> end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             !> simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             !> end
             !>  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(nn)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1         
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             !> simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             !> end
          enddo
       enddo
       nn=nn+1
    enddo
    enddo
    phi_ij = phi
    end subroutine Sweep2D
