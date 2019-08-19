   if (BC1 == 'Reflective Right') then
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    ! vacuum in the left side
    do k = 1,ng
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          phix(1,j,n) = 0.0d0  ! pour mu>0
          do i = 1,totNFM_x         
             phiy(i,1,n)     = 0.0d0 ! pour eta>0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n)   = 2.0d0*phi_ijn(i,j,n,k) - phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j+1,n)   = 2.0d0*phi_ijn(i,j,n,k) - phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do j = totNFM_y,1,-1
          if (BC1 == 'Reflective Right') then
          phix(totNFM_x+1,j,n) = phix(totNFM_x+1,j,numord+1-n) ! pour mu<0   
          else
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0   
          endif
          do i = totNFM_x,1,-1   
             phiy(i,totNFM_y+1,n) = 0.0D0 ! pour eta < 0     
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1 
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1
          phix(1,j,n) = 0.0D0             ! pour mu>0
          do i = 1,totNFM_x         
             phiy(i,totNFM_y+1,n) = 0.0D0 ! pour eta<0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j,n)   = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y
          if (BC1 == 'Reflective Right') then
          phix(totNFM_x+1,j,n) = phix(totNFM_x+1,j,numord+1-n) ! pour mu<0 
          else
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0
          endif 
          do i = totNFM_x,1,-1  
             phiy(i,1,n) = 0.0D0       ! pour eta>0    
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) =  2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    elseif (BC4 == 'Reflective Top') then
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    ! vacuum in the left side
    do k = 1,ng
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          phix(1,j,n) = 0.0d0  ! pour mu>0
          do i = 1,totNFM_x         
             phiy(i,1,n)     = 0.0d0 ! pour eta>0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n)   = 2.0d0*phi_ijn(i,j,n,k) - phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j+1,n)   = 2.0d0*phi_ijn(i,j,n,k) - phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do j = totNFM_y,1,-1
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0   
          do i = totNFM_x,1,-1   
             phiy(i,totNFM_y+1,n) = phiy(i,totNFM_y+1,numord+1-n) ! pour eta < 0     
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1 
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0
          do i = totNFM_x,1,-1  
             phiy(i,1,n) = 0.0D0       ! pour eta>0    
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) =  2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1
          phix(1,j,n) = 0.0D0             ! pour mu>0
          do i = 1,totNFM_x       
             phiy(i,totNFM_y+1,n) = phiy(i,totNFM_y+1,numord+1-n) ! pour eta<0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j,n)   = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    elseif (BC2 == 'Reflective Left') then
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do j = totNFM_y,1,-1
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0   
          do i = totNFM_x,1,-1   
             phiy(i,totNFM_y+1,n) = 0.0d0 ! pour eta < 0     
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1 
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    ! vacuum in the left side
    do k = 1,ng
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          phix(1,j,n) =  phix(1,j,numord+1-n)  ! pour mu>0
          
          do i = 1,totNFM_x         
             phiy(i,1,n)     = 0.0d0 ! pour eta>0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n)   = 2.0d0*phi_ijn(i,j,n,k) - phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j+1,n)   = 2.0d0*phi_ijn(i,j,n,k) - phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0
          do i = totNFM_x,1,-1  
             phiy(i,1,n) = 0.0D0       ! pour eta>0    
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) =  2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1
          phix(1,j,n) =  phix(1,j,numord+1-n)             ! pour mu>0
          do i = 1,totNFM_x       
             phiy(i,totNFM_y+1,n) = 0.0d0 ! pour eta<0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j,n)   = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    elseif (BC3 == 'Reflective Bottom') then
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do j = totNFM_y,1,-1
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0   
          do i = totNFM_x,1,-1   
             phiy(i,totNFM_y+1,n) = 0.0d0 ! pour eta < 0     
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1 
             phiy(i,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    ! vacuum in the left side
    do k = 1,ng
    do n = 3*numord/4+1,numord
       do j = 1,totNFM_y
          phix(1,j,n) =  0.0d0
          do i = 1,totNFM_x         
             phiy(i,1,n)     = phiy(i,1,numord+1-n) ! pour eta>0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n)   = 2.0d0*phi_ijn(i,j,n,k) - phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j+1,n)   = 2.0d0*phi_ijn(i,j,n,k) - phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/4+1,numord/2
       do j = totNFM_y,1,-1
          phix(1,j,n) =  0.0d0             ! pour mu>0
          do i = 1,totNFM_x       
             phiy(i,totNFM_y+1,n) = 0.0d0 ! pour eta<0
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i,j,n) + ay*phiy(i,j+1,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i+1,j,n) = 2.0D0*phi_ijn(i,j,n,k)-phix(i,j,n)
             ! simple negative flux fixup         
             if (phix(i+1,j,n)<0) phix(i+1,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = 1,totNFM_x
             phiy(i,j,n)   = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j+1,n)
             ! simple negative flux fixup         
             if (phiy(i,j,n)<0) phiy(i,j,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/2+1,3*numord/4
       do j = 1,totNFM_y
          phix(totNFM_x+1,j,n) = 0.0d0 ! pour mu<0
          do i = totNFM_x,1,-1  
             phiy(i,1,n) = phiy(i,1,numord+1-n)     ! pour eta>0    
             ax = 2.0d0*abs(mu(n))/dx(i)
             ay = 2.0d0*abs(eta(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*phix(i+1,j,n) + ay*phiy(i,j,n)
             deno = SigT(fmmid_xy(i,j),k) + ax + ay
             phi_ijn(i,j,n,k) = nume/deno
             ! simple negative flux fixup         
             if (phi_ijn(i,j,n,k)<0) phi_ijn(i,j,n,k) = 0.0d0
             ! end
             phix(i,j,n) =  2.0D0*phi_ijn(i,j,n,k)-phix(i+1,j,n)
             ! simple negative flux fixup         
             if (phix(i,j,n)<0) phix(i,j,n) = 0.0d0
             ! end
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*w(n)*phi_ijn(i,j,n,k)
          enddo
          do i = totNFM_x,1,-1
             phiy(i,j+1,n) = 2.0D0*phi_ijn(i,j,n,k)-phiy(i,j,n)
             ! simple negative flux fixup         
             if (phiy(i,j+1,n)<0) phiy(i,j+1,n) = 0.0d0
             ! end
          enddo
       enddo
    enddo
    enddo
    endif


!BL7
    subroutine Sweep(totNFM_x,totNFM_y,ngauss,Nmat,ng,Q_nij,dx,dy,SigT,&
                     mup,etap,pwi,fmmid_xy,phi_ij)
    implicit none
    integer(kind=4), intent(in) :: ngauss,Nmat,ng,totNFM_x,totNFM_y
    integer(kind=4), dimension(totNFM_x,totNFM_y), intent(in) :: fmmid_xy
    real(kind=8), dimension(totNFM_x,totNFM_y,ngauss*(ngauss+2)/2,ng), intent(in) :: Q_nij 
    real(kind=8), dimension(ngauss*(ngauss+2)/2), intent(in) :: mup, etap, pwi
    real(kind=8), dimension(totNFM_x), intent(in) :: dx
    real(kind=8), dimension(totNFM_y), intent(in) :: dy
    real(kind=8), dimension(Nmat,ng), intent(in) :: SigT
    real(kind=8), dimension(totNFM_x,totNFM_y,ng), intent(out) :: phi_ij
    real(kind=8), dimension(totNFM_y,ngauss*(ngauss+2)/2,ng) :: LeftBoundary, RightBoundary 
    real(kind=8), dimension(totNFM_x,ngauss*(ngauss+2)/2,ng) :: BottomBoundary, TopBoundary
    real(kind=8), dimension(ngauss*(ngauss+2)/2,ng) :: psiH, psiC
    real(kind=8), dimension(totNFM_y,ngauss*(ngauss+2)/2,ng) :: psiV
    real(kind=8), dimension(totNFM_x,totNFM_y,ng) :: phi
    real(kind=8) :: ax,ay,nume,deno
    integer(kind=4) :: i,j,n,k,numord
    LeftBoundary = 0.0D0; RightBoundary = 0.0D0
    BottomBoundary = 0.0D0; TopBoundary = 0.0D0
    psiH=0.0d0;psiV=0.0d0;psiC=0.0d0;phi=0.0d0
    numord = ngauss*(ngauss+2)/2
    !       ----------------------------------------------------------------- 
    !        mu > 0, eta > 0,   LEFT-TO-RIGHT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 3*numord/4+1,numord
       do i = 2,totNFM_x+1
          psiH(n,k) = BottomBoundary(i,n,k)
          do j = 2,totNFM_y+1
             ax = 2.0d0*abs(mup(n))/dx(i-1)
             ay = 2.0d0*abs(etap(n))/dy(j-1)
             nume = Q_nij(i-1,j-1,n,k) + ax*psiV(j-1,n,k) + ay*psiH(n,k)
             deno = ax + ay + SigT(fmmid_xy(i-1,j-1),k)
             psiC(n,k) = nume/deno
             if (psiC(n,k)<0) then
                 psiC(n,k)=0
             endif
             psiV(j-1,n,k) = 2.0D0*psiC(n,k)-psiV(j-1,n,k)
             psiH(n,k) = 2.0D0*psiC(n,k)-psiH(n,k)
             !  update the scalar flux
             phi(i-1,j-1,k) = phi(i-1,j-1,k) + 0.25D0*pwi(n)*psiC(n,k)
          enddo
       enddo
    enddo
    enddo  
    !       -----------------------------------------------------------------
    !        mu < 0, eta > 0,   RIGHT-TO-LEFT, BOTTOM-TO-TOP
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = numord/2+1,3*numord/4
       do i = totNFM_x,1,-1
          psiH(n,k) = BottomBoundary(i,n,k)
          do j = 2,totNFM_y+1
             ax = 2.0d0*abs(mup(n))/dx(i)
             ay = 2.0d0*abs(etap(n))/dy(j-1)
             nume = Q_nij(i,j-1,n,k) + ax*psiV(j-1,n,k) + ay*psiH(n,k)
             deno = ax + ay + SigT(fmmid_xy(i,j-1),k)
             psiC(n,k) = nume/deno
             if (psiC(n,k)<0) then
                 psiC(n,k)=0
             endif
             psiV(j-1,n,k) = 2.0D0*psiC(n,k)-psiV(j-1,n,k)
             psiH(n,k) = 2.0D0*psiC(n,k)-psiH(n,k)
             !  update the scalar flux
             phi(i,j-1,k) = phi(i,j-1,k) + 0.25D0*pwi(n)*psiC(n,k)
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu > 0, eta < 0,   LEFT-TO-RIGHT, TOP-TO-BOTTOM
    !       ----------------------------------------------------------------- 
    do k = 1,ng
    do n = numord/4+1,numord/2
       do i = 2,totNFM_x+1
          psiH(n,k) = TopBoundary(i,n,k)
          do j = totNFM_y,1,-1
             ax = 2.0d0*abs(mup(n))/dx(i-1)
             ay = 2.0d0*abs(etap(n))/dy(j)
             nume = Q_nij(i-1,j,n,k) + ax*psiV(j,n,k) + ay*psiH(n,k)
             deno = ax + ay + SigT(fmmid_xy(i-1,j),k)
             psiC(n,k) = nume/deno
             if (psiC(n,k)<0) then
                 psiC(n,k)=0
             endif
             psiV(j,n,k) = 2.0D0*psiC(n,k)-psiV(j,n,k)
             psiH(n,k) = 2.0D0*psiC(n,k)-psiH(n,k)
             !  update the scalar flux
             phi(i-1,j,k) = phi(i-1,j,k) + 0.25D0*pwi(n)*psiC(n,k)
          enddo
       enddo
    enddo
    enddo
    !       -----------------------------------------------------------------        
    !        mu < 0, eta < 0,   RIGHT-TO-LEFT, TOP-TO-BOTTOM
    !       -----------------------------------------------------------------
    do k = 1,ng
    do n = 1,numord/4
       do i = totNFM_x,1,-1
          psiH(n,k) = TopBoundary(i,n,k)
          do j = totNFM_y,1,-1
             ax = 2.0d0*abs(mup(n))/dx(i)
             ay = 2.0d0*abs(etap(n))/dy(j)
             nume = Q_nij(i,j,n,k) + ax*psiV(j,n,k) + ay*psiH(n,k)
             deno = ax + ay + SigT(fmmid_xy(i,j),k)
             psiC(n,k) = nume/deno
             if (psiC(n,k)<0) then
                 psiC(n,k)=0
             endif
             psiV(j,n,k) = 2.0D0*psiC(n,k)-psiV(j,n,k)
             psiH(n,k) = 2.0D0*psiC(n,k)-psiH(n,k)
             !  update the scalar flux
             phi(i,j,k) = phi(i,j,k) + 0.25D0*pwi(n)*psiC(n,k)
          enddo
       enddo
    enddo
    enddo
    phi_ij = phi
    end subroutine

!BL3
    subroutine LevelSymTable(mu,eta,w,N)
       implicit none
       integer(kind=4), intent(in) :: N
       real(kind=8), dimension(N*(N+2)/2), intent(out) :: mu,eta,w
       integer(kind=4) :: i,k1
       CHARACTER(10) :: ligne
       real(kind=8), parameter :: PI = 3.141592654
       real(kind=8) :: w1,w2
       open (10,file='app/data/ordinate.h') 
       if (N==2) then
       do
          ligne = ""
          read(10,"(a)") ligne 
          k1 = index( ligne , "2 " )
          if ( k1 /=0) exit
       enddo
       do i =1,N*(N+2)/2
          read(10,*) mu(i),eta(i)
       enddo
       w = 4.0D0*PI
       elseif (N==4) then
       do
          ligne = ""
          read(10,"(a)") ligne 
          k1 = index( ligne , "4 " )
          if ( k1 /=0) exit
       enddo
       do i =1,N*(N+2)/2
          read(10,*) mu(i),eta(i)
       enddo
       w = 1.0D0/3.0D0*PI*4
       elseif (N==6) then
       do
         ligne = ""
         read(10,"(a)") ligne 
         k1 = index( ligne , "6 " )
         if ( k1 /=0) exit
       enddo
       do i =1,N*(N+2)/2
       read(10,*) mu(i),eta(i)
       enddo
       w1 = 0.1761263; w2 = 0.1572071
       w(1:12)  = (/ w1, w2, w2, w1, w2, w1, w1, w2, w1, w1, w2, w1/)
       w(13:24) = w(12:1:-1)
       w=w*PI*4
       endif
    end subroutine LevelSymTable
