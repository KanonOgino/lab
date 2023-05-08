c**********************************************************************
c*    Check centroids H2 and H2O craster degree                       *
c*    Editor: K.Ogino 3/11/11                                         *
c**********************************************************************
      program naname_obz
c**********************************************************************

      implicit none
      integer :: i,j,k,l,n,m,nbead,natom,nstep
      integer :: x=1,y=2,z=3
      character (len=13) trj_name,outputfile
      character (len=2), allocatable :: atom(:)
      integer, allocatable :: atom_bead(:), step(:)
      double precision, allocatable :: atom_cart(:,:)
      double precision, allocatable :: h(:,:,:,:)
      double precision :: cent(1:3),cent_atom(1:3)
      double precision :: centatom(1:3)

c**********************************************************************
  
c**********************************************************************
c                          input file
c**********************************************************************

      namelist /input/ trj_name,outputfile,nstep,nbead,natom
      read(5,input)
c      print *, trj_name
c**********************************************************************
c**********************************************************************
c                          allocate
c**********************************************************************

      allocate(atom(natom))
      allocate(atom_cart(natom,3))
      allocate(atom_bead(nstep))
      allocate(step(nstep))
      allocate(h(nstep,nbead,natom,3))


c**********************************************************************

c**********************************************************************
c                          main
c**********************************************************************
      open(10, file=trim(trj_name), status='old')
      open(11, file=trim(outputfile), status='replace')

      !print *, trj_name
      do i = 1, nstep
        read(10,*,end=100) atom_bead(i)
        read(10,*) step(i)
        do j = 1, nbead
          do k = 1, natom
            read(10,*,end=100) atom(k),(atom_cart(k,l),l=1,3)
              !print *, atom(k)
              do l = 1,3
                h(i,j,k,l) = atom_cart(k,l)
              enddo
          enddo
        enddo
      enddo


      print *, 'half'

100     nstep = i-1 
        print *, nstep
        do i =  1, nstep
        !print *, nstep
        write(11,*) natom
        write(11,*) step(i)

        do k = 1, natom
          cent = 0.0d0
          cent_atom = 0.0d0
          do j = 1, nbead
            cent(x) = cent(x) + h(i,j,k,x)
          enddo

          do j = 1, nbead
            cent(y) = cent(y) + h(i,j,k,y)
          enddo

          do j = 1, nbead
            cent(z) = cent(z) + h(i,j,k,z)
          enddo

          cent_atom = cent/dble(nbead)

          write(11,*) atom(k), cent_atom

        enddo

!        cent_atom = cent_atom + cent/dble(nbead)

      enddo


      print *,'finish'

      close(10)
      close(11)
      deallocate(atom)
      deallocate(atom_cart)
      deallocate(atom_bead)
      deallocate(step)
      deallocate(h)
      end program naname_obz