module cliffs_tolkova_mod
    !!
    !!
    !! This module contains edited excerpts from the Cliffs solver. Cliffs was developed by 
    !! Elena Tolkova. It is based on the MOST solver but has an alternate treatment of 
    !! wetting and drying.
    !!
    !! The Cliffs code is available on github (https://github.com/Delta-function/cliffs-src)
    !! under a FreeBSD licence. The code here is a modified version of that. 
    !!
    !! Below we reproduce the copyright information of cliffs_main.f as required by the FreeBSD licence.
    !!
    !!

    !!! From cliffs_main.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! https://github.com/Delta-function/cliffs-src/blob/master/cliffs_main.f
    !! 
    !!  Cliffs is an open-source model for tsunami propagation and inundation computations       
    !!  in the framework of the shallow-water equations. Cliffs implements:                      
    !!  (1) VTCS-2 finite-difference scheme in an open ocean and on the open boundary            
    !! 	 as in (Titov and Synolakis, 1995, 1998);                                               
    !!  (2) dimensional splitting as in (Titov and Gonzalez, 1997; Titov and Synolakis, 1998);   
    !!  (3) reflection and inundation computations as in (Tolkova, 2014);                        
    !!  (4) data flows similar to curvilinear MOST and MOST-4 (Burwell and Tolkova, 2008).       
    !!                               REFERENCE:                                                  
    !!  E. Tolkova. Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting 
    !!  Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314                     
    !!                                                                                           
    !!  FEATURES: Cartezian or Geophysical (lon/lat) coordinates; 2D or 1D domains;              
    !!           grid nesting with one-way coupling; initial conditions or/and boundary forcing;
    !!           Open MP; NetCDF format of I/O                                                  
    !!
    !!           Copyright (c) 2014, Elena Tolkova                                              
    !!           under the terms of FreeBSD License                                             
    !!                                                                                          
    !!  Redistribution and use in source and binary forms, with or without                       
    !!  modification, are permitted provided that the following conditions are met:              
    !!                                                                                          
    !!   1. Redistributions of source code must retain the above copyright notice, this          
    !!     list of conditions and the following disclaimer.                                     
    !!   2. Redistributions in binary form must reproduce the above copyright notice,            
    !!     this list of conditions and the following disclaimer in the documentation            
    !!     and/or other materials provided with the distribution.                               
    !!                                                                                          
    !!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND          
    !!  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED            
    !!  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE                   
    !!  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR          
    !!  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES           
    !!  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;             
    !!  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND              
    !!  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT               
    !!  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS            
    !!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                             
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! From cliffs.f !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! https://github.com/Delta-function/cliffs-src/blob/master/cliffs.f
    !!
    !! Cliffs SOLVER to propagate tsunami wave through one time step in 1D with:                !
    !! (1) VTCS-2 finite-difference scheme in an open ocean and on the open boundary            !
    !!	 as in (Titov and Synolakis, 1995, 1998)                                             !
    !! (2) reflection and inundation computations as in (Tolkova, 2014)                         !
    !!                               REFERENCES:                                                !   
    !! V. Titov and C. Synolakis. Modeling of Breaking and Nonbreaking Long-Wave Evolution and  !
    !! Runup Using VTCS-2. J. of Waterway, Port, Coastal, and Ocean Eng. Vol. 121, No 6 (1995), ! 
    !! pp. 308-316.                                                                             !
    !! E. Tolkova. Land-Water Boundary Treatment for a Tsunami Model With Dimensional Splitting.!
    !! Pure and Applied Geophysics, Vol. 171, Issue 9 (2014), pp. 2289-2314                     !
    !!                                                                                          !
    !!     Copyright (C) 2014, Elena Tolkova                                                    !
    !!     For conditions of distribution and use, see copyright notice in cliffs_main.f        !
    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! SWALS parameters
    use global_mod, only : dp, ip, gravity
    use logging_mod, only: log_output_unit

    implicit none

    ! Parameters required by Cliffs
    real(dp), parameter :: grav = gravity, sqr2 = sqrt(2.0_dp)

    contains

    subroutine cliffs(dt, iu, irw, n, cliffs_minimum_allowed_depth, &
                      elev, cel, xvel, yvel, s1, s2, zeta, manning_squared, &
                      linear_friction_coeff)
        !! This is the Cliffs 1D time-stepping routine. 
        !! It was originally taken from the cliffs source-code by Elena Tolkova, see comments above for licence etc.
        !! Compared with the cliffs original version, the code below has various interface edits to cleanly integrate into SWALS.

        ! GD comment -- arguments are:
        !    dt -- timestep
        !    iu -- integer, has values 1 or 2 -- referring to computations in x(lon) or y(lat) slice
        !    irw -- the y or x index of the grid that this slice represents
        !    n -- the number of x or y values in this slice.
        !    cliffs_minimum_allowed_depth -- the 'ground' parameter in cliffs, which is a wet/dry depth threshold.
        !                                    This value needs to be tuned for the problem.
        ! GD Comments
        !
        ! elev(nXn,nYn) -- used in place of 'dep' in cliffs.f, which is the "undisturbed water depth" or "depth below MSL" -- 
        !                  same as my "-domain%U(,,ELV)". So below the original cliffs routine is modified to 
        !                  use '-elev' instead of 'dep'
        ! cel(nXn,nYn) -- sqrt(g * real_water_depth). Note the code sometimes uses 'h' to represent this 
        ! xvel(nXn,nYn) -- x velocity
        ! yvel(nXn,nYn) -- y velocity
        ! s1 -- x-direction cell size -- for SWALS, I reduced it to a constant (for each slice)
        ! s2 -- y-direction cell size -- for SWALS I reduced it to a constant (for each slice)
        ! zeta(nYn) -- spherical coordinate scale factor. In CLIFFS this is 0.125/Earth_Radius * tan(pi/180 * lat)
        ! manning_squared -- like 'crough' in the original cliffs.f source -- except for the version below, by default we use
        !                    the standard manning formula (despite the expensive power-law computation)
        ! linear_friction_coeff -- drag term, which if written in terms of the depth-integrated-velocity equations, would be
        !                          like subtracting "linear_friction_coef * UH (or VH)" from the right-hand-side of the UH/VH 
        !                          evolution equations. This form of drag is used in Fine et al (2012) and Kulikov et al (2014)
        real(dp), intent(in):: dt, cliffs_minimum_allowed_depth
        integer(ip), intent(in) :: iu,irw,n
        real(dp), intent(in) :: elev(:,:), s1, s2, zeta(:), manning_squared(:,:) !, edge1(:,:), edge2(:,:)
        real(dp), intent(inout) :: cel(:,:), xvel(:,:), yvel(:,:)
        real(dp), intent(in) :: linear_friction_coeff

        real(dp) :: h(n),pp(n),qq(n),vv(n),dx(n),scl(n)
        real(dp) :: depth(n),Pinv(n),Qinv(n),v(n),u(n), crough(n)
        real(dp) :: pqvdL(4,99), pqvdR(4,99)
        real(dp) :: uj,vj,pj,qj,dpj,dpm,dtx,hj
        real(dp) :: ej,e2,e1,Qpm,Qpj,Qjm,Qsum
        real(dp) :: v1,v2,d1,d2,p1,p2,q1,q2      
        real(dp) :: cc,sphr,fm,fp,flood     
        integer(ip) :: i,j,i1,i2,kseg,k,k1,k2
        integer(ip) :: lghost(n),rghost(n)
        logical :: land(n),water,newland(n)
        real(dp) :: ground, celmin

        ! Cliffs parameters
        ground = cliffs_minimum_allowed_depth
        celmin = sqrt(gravity * cliffs_minimum_allowed_depth)


        ! GD -- This is packing a 1D slice of the 2D grid into a 1D array (the algorithm is based on dimension splitting, so is
        !       essentially 1D)
        if(iu.eq.1) then
            ! GD -- iu == 1, Computations in lon (x) direction
             do i = 1,n
                depth(i)= -elev(i,irw) ! GD Comment -- here we use negative elevation in place of 'dep' in CLIFFS
                h(i) = cel(i,irw)
                u(i) = xvel(i,irw)
                v(i) = yvel(i,irw)
                dx(i) = s1
                scl(i) = 0
                crough(i) = manning_squared(i, irw)
            end do 
        else
            ! GD -- iu == 2, Computations in lat (y) direction
            do i = 1,n
                depth(i)= -elev(irw,i) ! GD Comment -- here we use negative elevation in place of 'dep' in CLIFFS
                h(i) = cel(irw,i)
                u(i) = yvel(irw,i)
                v(i) = xvel(irw,i)
                dx(i) = s2
                scl(i) = zeta(i) 
                crough(i) = manning_squared(irw, i)
           end do
        endif


        !
        ! GD -- note below here, everything is 1D, until the very end where arrays are repacked.
        !

        lghost(1:n)=0
        rghost(1:n)=0
        ! Find dry areas
        land=.false.
        do j=1,n 
            if(h(j).le.celmin) land(j)=.true.
        end do
        newland=land
        
        ! GD -- figure out if the left-most-cell is land or water 
        water=.false.
        kseg=0
        if(.not.land(1)) then
            water=.true.
            kseg=1
            lghost(1)=1
        endif
        k1=kseg+1 ! first segment with left shore
          
        ! GD -- some wetting and drying stuff 
        do i= 2,n
            if(water) then
                ! Here segment (i-1) is wet
                if(land(i)) then
                    ! Here segment (i) is dry

                    ! GD Comment -- notice h**2 / grav = (grav * real_water_depth / grav ) = real_water_depth.
                    !               Whereas 'depth' here refers to the "undisturbed water depth"
                    !               So "(h(i-1)**2)/grav-depth(i-1)" is like stage(i-1),
                    !               flood = stage(i-1) + undisturbed_water_depth(i)
                    flood=(h(i-1)**2)/grav-depth(i-1)+depth(i)
                    if((flood.gt.ground).and.(.not.land(i-1)).and.(u(i-1).gt.0)) then ! expand
                        ! GD Comment -- here, the dry cell has u/v/depth extrapolated from its left neighbour,
                        !               because 'flood' was deep enough to suggest flooding occurs
                        u(i)=u(i-1)
                        v(i)=v(i-1)
                        h(i)=celmin
                        newland(i)=.false.
                    else
                        ! GD Comment -- here, segment (i) remains dry
                        water=.false. ! close wet segment
                        rghost(kseg)=i
                    endif
                endif
            else
                ! Here segment(i-1) is dry
                if(.not.land(i)) then ! start wet segment
                    ! Here segment (i) is wet
                    water = .true.
                    kseg = kseg+1
                    ! GD Comment -- so 'flood' = stage(i) + undisturbed_water_depth(i-1)
                    !               See comment above for details
                    flood = (h(i)**2)/grav-depth(i)+depth(i-1) 
                    if((flood.gt.ground).and.(u(i).lt.0)) then ! expand
                        !! Here the initially dry 'i-1' cell has u/v/depth extrapolated because it
                        !! seems flooding should occur, based on 'flood'
                        lghost(kseg)=max(i-2,1)
                        newland(i-1)=.false.
                        u(i-1)=u(i)
                        v(i-1)=v(i)
                        h(i-1)=celmin
                    else
                        !! Here segment (i-1) remains dry
                        lghost(kseg)=i-1
                    end if

                    if(kseg.gt.1) then
                        if((lghost(kseg).lt.rghost(kseg-1))) then ! glue to previous wet segment 
                            lghost(kseg)=0
                            kseg=kseg-1
                            rghost(kseg)=0
                        end if
                    end if
                endif
            endif
        end do
          
        if(kseg.eq.0) return !goto 99
        if (rghost(kseg).eq.0) rghost(kseg)=n
        k2=kseg-1
        if(newland(rghost(kseg))) k2=kseg
          
        where(newland)
            h=0
            u=0
            v=0
        end where
        vv=v 

        !   compute Riemann invarients
        Pinv=0
        Qinv=0
        !   - in open sea	
        do k=1,kseg
            i1=lghost(k)+1
            i2=rghost(k)-1
            do i = i1,i2
                Pinv(i)=u(i)+2*h(i)
                Qinv(i)=u(i)-2*h(i)
            end do
        enddo
        !    - on left shoreline in ghost node
        if (k1.le.kseg) then 
            do k=k1,kseg
                i=lghost(k)
                j=i+1
                pqvdL(1:4,k)=(/-Qinv(j), -Pinv(j), v(j), depth(j)/)         
            enddo
        endif 
        !    - on right shoreline in ghost node
        if (k2.ge.1) then
            do k=1,k2
                i=rghost(k)
                j=i-1
                pqvdR(1:4,k)=(/-Qinv(j), -Pinv(j), v(j), depth(j)/)         
            end do
        endif 
        !    - on left wet edge
        if(.not.newland(1)) then
            Pinv(1)=u(1)+2*h(1)
            Qinv(1)=u(1)-2*h(1)
            pqvdL(1:4,1)=(/Pinv(1), Qinv(1), v(1), depth(1)/)
        end if
        !    - on right wet edge
        if(.not.newland(n)) then
            Pinv(n)=u(n)+2*h(n)
            Qinv(n)=u(n)-2*h(n)
            pqvdR(1:4,kseg)=(/Pinv(n), Qinv(n), v(n), depth(n)/)
        end if
        qq=Qinv
        pp=Pinv
                      
        ! time-step integration      
        do k=1,kseg
     
            i1=lghost(k)+1
            i2=rghost(k)-1
            if(i1.gt.i2) cycle !goto 9
            
            do i=i1,i2

                j=i-1
                if(i==i1) then
                    p1=pqvdL(1,k)
                    q1=pqvdL(2,k)
                    v1=pqvdL(3,k)
                    d1=pqvdL(4,k)
                else
                    p1=Pinv(j)
                    q1=Qinv(j)
                    v1=v(j)
                    d1=depth(j)
                end if
                j=i+1;
                if(i==i2) then
                    p2=pqvdR(1,k)
                    q2=pqvdR(2,k)
                    v2=pqvdR(3,k)
                    d2=pqvdR(4,k)
                else
                    p2=Pinv(j)
                    q2=Qinv(j)
                    v2=v(j)
                    d2=depth(j)
                end if
                uj=u(i)
                vj=v(i)
                ! GD Comment -- hj = real_water_depth(i)
                hj=h(i)**2/grav
                !  Almost Manning - Preferred friction model     
                !cc=grav*crough(i)*uj*sqrt((uj**2+vj**2)/hj)/hj  
                !   Drag force 
                !  		cc=grav*crough(i)*uj*sqrt(uj**2+vj**2)/hj
                !  Manning friction - SLOW -- but for SWALS we use this by default
                cc = grav*crough(i)*uj*sqrt(uj**2+vj**2)/hj**(4.0_dp/3.0_dp) + &
                    ! Here we append a linear drag model like Fine et al., (2012), Kulikov et al., (2014).
                    linear_friction_coeff * uj
                
                pj=Pinv(i)
                qj=Qinv(i)
                sphr=(pj-qj)*(pj+qj)*scl(i)
                
                dpj=grav*(d2-depth(i))
                dpm=grav*(d2-d1)        
                fm=dx(i-1)
                fp=dx(i)
                dtx=dt/(fm+fp)
                
                ej=3*pj+qj
                e2=3*p2+q2
                e1=3*p1+q1
                Qpm=0.125*(e2+e1)*(p2-p1)-dpm 
                Qpj=0.125*(e2+ej)*(p2-pj)-dpj 
                Qjm=0.125*(ej+e1)*(pj-p1)-dpm+dpj
                Qsum=Qpm+0.25*ej*dt*(Qjm/fm-Qpj/fp)
                pp(i)=pj-dtx*Qsum-dt*(cc-sphr)
                
                ej=3*qj+pj
                e2=3*q2+p2
                e1=3*q1+p1
                Qpm=0.125*(e2+e1)*(q2-q1)-dpm
                Qpj=0.125*(e2+ej)*(q2-qj)-dpj
                Qjm=0.125*(ej+e1)*(qj-q1)-dpm+dpj
                Qsum=Qpm+0.25*ej*dt*(Qjm/fm-Qpj/fp)
                qq(i)=qj-dtx*Qsum-dt*(cc+sphr)
                
                Qsum=(v2-vj)/fp-(vj-v1)/fm
                vv(i)=vj-uj*dtx*(v2-v1-Qsum*(uj*dt+fp-fm))
            end do
        end do

        do k=1,kseg
            i1=lghost(k)+1
            i2=rghost(k)-1
            do i=i1,i2
                h(i)=(pp(i)-qq(i))/4
                u(i)=(pp(i)+qq(i))/2
            end do
        end do
          
        !   check water surface angle on left shoreline
        if (k1.le.kseg) then 
            do k=k1,kseg
                i=lghost(k)+1
                j=i+1
                if(land(i)) then ! if newly included node
                    if(sqr2*h(i).gt.h(j)) then
                        h(i)=h(j)/sqr2
                        u(i)=u(j)/2
                        vv(i)=vv(j)/2
                    endif
                endif
            enddo
        endif 
        !   check water surface angle on right shoreline
        if (k2.ge.1) then
            do k=1,k2
                i=rghost(k)-1
                j=i-1
                if(land(i)) then ! if newly included node
                    if(sqr2*h(i).gt.h(j)) then
                        h(i)=h(j)/sqr2
                        u(i)=u(j)/2
                        vv(i)=vv(j)/2
                    endif
                endif
            enddo
        endif
         
        !! For SWALS we can skip these boundary conditions, and rely on the regular boundary treatment. 
        !!
        !! Integration in edges
        !flood=edge1(irw,3)+depth(1)
        !if((.not.newland(1)).and.(.not.newland(2)) .and.(flood.gt.ground)) then
        !    pj=Pinv(1)
        !    qj=Qinv(1)
        !    p2=Pinv(2)
        !    q2=Qinv(2)
        !    sphr=(pj-qj)*(pj+qj)*scl(1)
        !    Qpj=(3*(q2+qj)+p2+pj)*(q2-qj)/8-grav*(depth(2)-depth(1))
        !    qq(1)=qj-dt*(Qpj/dx(1)+sphr)
        !    pp(1)=edge1(irw,iu)+2*sqrt(grav*flood)
        !    if(u(1).lt.0) then
        !        vv(1)=v(1)-dt*(v(2)-v(1))*u(1)/dx(1) 
        !    else
        !        vv(1)=edge1(irw,3-iu)
        !    end if
        !    h(1)=(pp(1)-qq(1))/4
        !    u(1)=(pp(1)+qq(1))/2
        !endif
        !flood=depth(n)+edge2(irw,3)
        !if((.not.newland(n-1)).and.(.not.newland(n)) .and.(flood.gt.ground)) then
        !    pj=Pinv(n)
        !    qj=Qinv(n)
        !    p1=Pinv(n-1)
        !    q1=Qinv(n-1)
        !    sphr=(pj-qj)*(pj+qj)*scl(n)
        !    Qjm=(3*(p1+pj)+q1+qj)*(pj-p1)/8-grav*(depth(n)-depth(n-1)) 
        !    pp(n)=pj-dt*(Qjm/dx(n-1)-sphr)   
        !    qq(n)=edge2(irw,iu)-2*sqrt(grav*flood)

        !    if(u(n).gt.0) then
        !       vv(n)=v(n)-dt*(v(n)-v(n-1))*u(n)/dx(n-1) 
        !    else
        !       vv(n)=edge2(irw,3-iu)
        !    end if
        !    h(n)=(pp(n)-qq(n))/4
        !    u(n)=(pp(n)+qq(n))/2
        !endif

        !
        ! GD comment -- here the arrays are re-packed.
        !
          
        if(iu.eq.1_ip) then
            do i = 1,n
                cel(i,irw)=h(i)
                xvel(i,irw)=u(i)
                yvel(i,irw)=vv(i)
            end do
        else
            do i=1,n
               cel(irw,i)=h(i)
               xvel(irw,i)=vv(i)
               yvel(irw,i)=u(i)
            end do
        endif        
          
        ! 99  	continue
    end subroutine cliffs


    subroutine setSSLim(nXn,nYn,bbb,alpha,depmin)
        !! The Cliffs bathymetry smoother. 
        !! This bathymetry smoothing routine is from the CLIFFS source code by Elena Tolkova. It also has superficial edits to cleanly
        !! integrate into SWALS. See comments at the top of this file for licence information etc.
        !! Cliffs tends to go unstable if the bathymetry is insufficiently smooth (see analysis by Tolkova in the CLIFFS manual), so
        !! this routine is very important.
    
        integer(ip), intent(in) :: nXn,nYn
        real(dp), intent(inout) :: bbb(nXn,nYn)
        real(dp), intent(in) :: alpha,depmin

        real(dp) :: c1,c2,alp2,dep
        real(dp), dimension(:), allocatable :: q,p
        logical, dimension(:), allocatable :: land
        integer(ip) nodes,n,i,j,itry,nm

        nm=max(nXn,nYn)
        allocate(q(nm),p(nm))
        allocate(land(nm))

        write(log_output_unit,*) 'Applying CLIFFS bathymetry smoother with alpha=', alpha
        write(log_output_unit,*) '    Node-wise depth adjustments in open water:'
        ! Repeat the smoother a fixed number of times
        do itry=1,12

            nodes=0
            ! Smooth in y direction
            do i=1,nXn
                q(1:nYn)=bbb(i,1:nYn)
                p(1:nYn)=q(1:nYn)
                land(1:nm)=.false.
                do j=1,nYn
                    if (q(j).lt.depmin) then
                        land(j)=.true.
                    else
                        q(j)=sqrt(q(j))
                    endif
                enddo
                n=nYn-1
                do j=2,n
                    if(.not.(land(j-1).or.land(j).or.land(j+1))) then
                        c1=q(j+1)/q(j)
                        c2=q(j-1)/q(j)
                        if ((c1.gt.alpha).or.(c2.gt.alpha)) then
                            q(j)=(q(j-1)+q(j+1))/2
                            p(j)=q(j)*q(j)
                            nodes=nodes+1
                        endif
                    endif
                enddo
                bbb(i,1:nYn)=p(1:nYn)
            enddo

            ! Smooth in x direction
            do i=1,nYn
                q(1:nXn)=bbb(1:nXn,i)
                p(1:nXn)=q(1:nXn)
                land(1:nm)=.false.
                do j=1,nXn
                    if (q(j).lt.depmin) then
                        land(j)=.true.
                    else
                        q(j)=sqrt(q(j))
                    endif
                enddo
                n=nXn-1
                do j=2,n
                    if(.not.(land(j-1).or.land(j).or.land(j+1))) then
                        c1=q(j+1)/q(j)
                        c2=q(j-1)/q(j)
                        if ((c1.gt.alpha).or.(c2.gt.alpha)) then
                            q(j)=(q(j-1)+q(j+1))/2
                            p(j)=q(j)*q(j)
                            nodes=nodes+1
                        endif
                    endif
                enddo
                bbb(1:nXn,i)=p(1:nXn)
            enddo
            write(log_output_unit,*) '    passage',itry,':',nodes,' nodes adjusted'
            if(nodes.eq.0) exit

        enddo  ! end itry loop

        if(nodes.ne.0) write(log_output_unit,*) '    ATTN: processing is not complete - too rough bathy or too small alpha'

        !
        ! A second-type of smoothing below here
        !
            
        write(log_output_unit,*) '    Adjustments at coastlines:'
        alp2=alpha**2
        do itry=1,12
        
            nodes=0
            ! Along y direction
            do i=1,nXn
                q(1:nYn)=bbb(i,1:nYn)
                land(1:nm)=.false.
                do j=1,nYn
                    if (q(j).lt.depmin) land(j)=.true.
                enddo
                n=nYn-1
                do j=2,n
                    if(.not.land(j)) then
                        if(land(j+1).and.(.not.land(j-1)) .or. (.not.land(j+1)).and.land(j-1)) then
                            if(land(j+1)) then
                                dep=q(j-1)
                            else
                                dep=q(j+1)
                            endif
                            if( dep/q(j).gt.alp2 ) then
                                q(j)=dep/alp2
                                nodes=nodes+1
                            endif
                        endif
                    endif
                enddo
                bbb(i,1:nYn)=q(1:nYn)
            enddo

            ! Along x direction
            do i=1,nYn
                q(1:nXn)=bbb(1:nXn,i)
                land(1:nm)=.false.
                do j=1,nXn
                    if (q(j).lt.depmin) land(j)=.true.
                enddo
                n=nXn-1
                do j=2,n
                    if(.not.land(j)) then
                        if(land(j+1).and.(.not.land(j-1)) .or.(.not.land(j+1)).and.land(j-1)) then
                            if(land(j+1)) then
                                dep=q(j-1)
                            else
                                dep=q(j+1)
                            endif
                            if( dep/q(j).gt.alp2 ) then
                                q(j)=dep/alp2
                                nodes=nodes+1
                            endif
                        endif
                    endif
                enddo
                bbb(1:nXn,i)=q(1:nXn)
            enddo
            write(log_output_unit,*) '    passage',itry,':',nodes,' nodes adjusted'
                    
            if(nodes.eq.0) exit
        enddo  ! end itry loop
        
        if(nodes.ne.0) write(log_output_unit,*) '    ATTN: processing is not complete - too rough bathy or too small alpha'

    end subroutine setSSLim

end module
