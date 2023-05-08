
c**********************************************************************
      program calc_ene_ice
c**********************************************************************

      implicit none
      integer :: i,j,k,l,natom,nbead
      character (len=8), allocatable :: atom(:)
      double precision, allocatable :: a_c(:,:,:,:),cn(:,:,:)
      double precision :: m(26),mtrs,sectrs,angbohr,dt,ene_h2,ene_ice
      double precision :: m_H = 1.00782504, m_O = 15.9949146
      double precision :: c_h(2,3),c_ice(2,3),v_v_h(3),v_v_ice(3),
     &mh_tot,mi_tot,norm_h,norm_ice

      natom = 26
      nbead = 96

c*****unit transform***************************************************
c unit transform      
c**********************************************************************
      mtrs = 1822.89d0 !atomic_weight to atomic_unit
      sectrs =  41.34138d0 !femtosecond to atomic_unit 
      angbohr = 1.889726d0 !angstrom to atomic_unit

      m_H = m_H * mtrs
      m_O = m_O * mtrs
      dt = 0.1d0 * sectrs
c**********************************************************************


      allocate(a_c(2,nbead,natom,3))
      allocate(atom(natom))
      allocate(cn(2,natom,3))
      open(10,file="trj.xyz",status='old')


c**********************************************************************
c read trj.xyz file
c**********************************************************************
      do i = 1,2
        do j = 1,nbead
          do k = 1,natom
            do l = 1,3
              a_c(i,j,k,l) = 0.0d0
            enddo
          enddo
        enddo
      enddo
      
      do i = 1,2
        read(10,*)
        read(10,*)
        do j = 1,nbead
          do k = 1,natom
            read(10,*) atom(k),  (a_c(i,j,k,l),l=1,3)
          enddo
        enddo
      enddo

c      print *, a_c

c**********************************************************************

c**********************************************************************
c difine m_tot
c**********************************************************************

      do i = 1, natom
        m(i) = 0.0d0
      enddo
      mh_tot = 0.0d0
      mi_tot = 0.0d0

      do i = 1,natom
        if(i .le. 18) then
          m(i) = m_H
        else
          m(i) = m_O
        endif
      enddo

      do i = 1,2
        mh_tot = mh_tot + m(i)
      enddo
      do i = 3,26
        mi_tot = mi_tot + m(i)
      enddo

c      print *, 'test',mh_tot,mi_tot 

c**********************************************************************

c**********************************************************************
c centroids of each beads (cn)
c**********************************************************************

      do i = 1,2
        do k = 1,natom
          do l = 1,3
            cn(i,k,l) = 0.0d0
          enddo
        enddo
      enddo

      do i = 1,2
        do j = 1,nbead
          do k = 1,natom
            do l = 1,3
             cn(i,k,l) = cn(i,k,l) + a_c(i,j,k,l) 
            enddo
          enddo
        enddo
      enddo

      do i = 1,2
        do j = 1,natom
          do k = 1,3
            cn(i,j,k) = cn(i,j,k)/dble(nbead)
          enddo
        enddo
      enddo

c      print *,'test2', cn

c**********************************************************************
 
c**********************************************************************
c centroids of each moleculas
c**********************************************************************

      do i = 1,2
        do j = 1,3
          c_h(i,j) = 0.0d0
          c_ice(i,j) = 0.0d0
        enddo
      enddo

      do i = 1,2
        do k = 1,3
          do j = 1,2
            c_h(i,k) = c_h(i,k) + (cn(i,j,k) * m(j))
          enddo
          do j = 3,26
            c_ice(i,k) = c_ice(i,k) + (cn(i,j,k) * m(j))
          enddo
        enddo
      enddo

      do i = 1,2
        do j = 1,3
          c_h(i,j) = c_h(i,j)/mh_tot
          c_ice(i,j) = c_ice(i,j)/mi_tot
        enddo
      enddo

c     print *, 'test3' , c_h, c_ice
c      print *, c_ice,mi_tot

c**********************************************************************

c**********************************************************************
c make vectors of velocity
c**********************************************************************

      do i =1,3
        v_v_h(i) = 0.0d0
        v_v_ice(i) = 0.0d0
      enddo


      do i = 1,3
        v_v_h(i) = c_h(2,i) - c_h(1,i) 
      enddo
      do i = 1,3
        v_v_ice(i) = c_ice(2,i) - c_ice(1,i)
      enddo

      do i = 1,3
        v_v_h(i) = v_v_h(i)/dt
      enddo
      do i = 1,3
        v_v_ice(i) = v_v_ice(i)/dt
      enddo

      do i = 1,3
        norm_h = norm_h + (v_v_h(i)*angbohr)**2
        norm_ice = norm_ice + (v_v_ice(i)*angbohr)**2
      enddo

      norm_h = sqrt(norm_h)
      norm_ice= sqrt(norm_ice)

c      print *, 'test', v_v_h, v_v_ice

c**********************************************************************
c calculation Energy
c**********************************************************************

      ene_h2 = 0.0d0
      ene_ice = 0.0d0

      do i = 1,3
        ene_h2 = ene_h2 + (0.5d0 * mh_tot * (v_v_h(i) * angbohr)**2)
      enddo

      do i = 1,3
        ene_ice = ene_ice +
     &(0.5d0 * mi_tot * (v_v_ice(i) * angbohr)**2)
      enddo

      ! print *, v_v_h
      ! print *, v_v_ice
      ! print *, mi_tot, mh_tot
c**********************************************************************

c**********************************************************************
c output
c**********************************************************************

      write(*,1000) ene_h2,ene_ice,norm_h,norm_ice
1000  format(e21.12,',',e21.12,',',e21.12,',',e21.12)
c**********************************************************************
      end
