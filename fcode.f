c *** This is a test program to drive the library functionality form the fortran
c *** side. It produces an all_susrf-to-all_surf2 dependence. When it runs, it
c *** should produce a constant value for the data transfered to the surface1
c *** elements, corresponding to the formula: val = sum_p (p*num_elem_surf2)
c *** IN 2017/05/01

      Program face_test
      use MPI
      Implicit none
      Integer :: im,jm,nel1,nel2,nno
      Integer,allocatable,dimension(:,:) :: icon,jcon
      Integer,allocatable,dimension(:) :: iacc
      Real*8, allocatable,dimension(:,:) :: xi,xj
      Integer icomw,irank,nproc,ier

      Integer ihandle, num_recv, num_send
      Integer,allocatable,dimension(:) :: irdis,ircnt, isdis,iscnt
      Integer,allocatable,dimension(:) :: irecv,isend
      Real*8, allocatable,dimension(:) :: area

      icomw = MPI_COMM_WORLD
      call MPI_INIT(ier)
      call MPI_COMM_SIZE(icomw, nproc, ier)
      call MPI_COMM_RANK(icomw, irank, ier)

      IF( irank.eq.0 ) PRINT*,'Executing from FORTRAN'

      im = 2
      jm = 3
      nel1 = im*jm
      nel2 = im*jm*2
      nno = (im+1)*(jm+1)

      allocate(icon(0:4,nel1), xi(3,nno))
      allocate(jcon(0:4,nel2), xj(3,nno))
      allocate(iacc(0:nproc-1))

c---- create two surfaces
      call make_faces( im,jm, icon, xi, jcon, xj )
c---- offset the node coefficients by process ID; identical nodes for both
      xi(1,:) = xi(1,:) + dble(irank*im)
      xj(1,:) = xj(1,:) + dble(irank*im)
c---- display surfaces with tecplot/paraview
      call dump_faces( irank, im,jm, icon, xi, jcon, xj )

c---- make an acceleration array; an array that can filter out overlaps by proc
      iacc(:) = nel2    ! positive values imply possible overlap with process

c---- invoke the face-matching library
      call incg_facematch_setOverlapFunction()

      call incg_Facematch_Init( ihandle, icomw,
     &                          nel1, icon, xi,
     &                          nel2, jcon, xj, iacc, ier )

      call incg_Facematch_GetSizes( ihandle, num_recv, num_send, ier )
      IF( irank.eq.0 ) PRINT*,'GOT SIZES', num_recv, num_send

      allocate( irdis(0:nproc-1), ircnt(0:nproc-1) )
      allocate( isdis(0:nproc-1), iscnt(0:nproc-1) )
      allocate( irecv(num_recv), isend(num_send) )
      allocate( area(num_recv) )
      irecv = -10     ! these numbers are dummy for trapping possible errors
      isend = -20     ! these numbers are dummy for trapping possible errors
      area = -9.9d0   ! these numbers are dummy for trapping possible errors

      call incg_Facematch_FillArrays( ihandle,
     &                                isend, irecv, area,
     &                                irdis, ircnt, isdis, iscnt, ier )

c *** At this point the arrays (on the Fortrna side) have been filled with data

      call incg_Facematch_Term( ihandle, ier )

c *** Here we perform a test of the communication functionality

      call match_test( icomw, irank, nproc, im,jm, icon, xi, jcon, xj,
     &                 num_recv, num_send, irecv, isend, area,
     &                 irdis, ircnt, isdis, iscnt )

      call MPI_FINALIZE(ier)

      deallocate( icon,jcon,iacc,xi,xj )
      deallocate( irdis, ircnt, isdis,iscnt, irecv,isend,area )

      End

      Subroutine make_faces( im,jm, icon, xi, jcon, xj )
      Implicit none
      Integer im,jm
      Integer,dimension(0:4,im*jm) :: icon
      Integer,dimension(0:4,im*jm*2) :: jcon
      Real*8, dimension(3,(im+1)*(jm+1)) :: xi,xj
      Integer i,j,n


      Do j = 1,jm+1
      Do i = 1,im+1
         n = (j-1)*(im+1) + i

         xi(1,n) = dble(i)
         xi(2,n) = dble(j)
         xi(3,n) = 0.0d0

         xj(1,n) = dble(i)
         xj(2,n) = dble(j)
         xj(3,n) = 0.1d0
      Enddo
      Enddo

      Do j = 1,jm
      Do i = 1,im
         n = (j-1)*im + i

         icon(0,n) = 4
         icon(1,n) = (j-1)*(im+1) + (i  )
         icon(2,n) = (j-1)*(im+1) + (i+1)
         icon(3,n) = (j  )*(im+1) + (i+1)
         icon(4,n) = (j  )*(im+1) + (i  )
      Enddo
      Enddo

      n = 0
      Do j = 1,jm
      Do i = 1,im
         n = n + 1

         jcon(0,n) = 3
         jcon(1,n) = (j-1)*(im+1) + (i  )
         jcon(2,n) = (j-1)*(im+1) + (i+1)
         jcon(3,n) = (j  )*(im+1) + (i  )

         n = n + 1

         jcon(0,n) = 3
         jcon(1,n) = (j  )*(im+1) + (i  )
         jcon(2,n) = (j-1)*(im+1) + (i+1)
         jcon(3,n) = (j  )*(im+1) + (i+1)
         jcon(4,n) = jcon(3,n)
      Enddo
      Enddo

      return
      end

      Subroutine dump_faces( irank, im,jm, icon, xi, jcon, xj )
      Implicit none
      Integer irank, im,jm
      Integer,dimension(0:4,im*jm) :: icon
      Integer,dimension(0:4,im*jm*2) :: jcon
      Real*8, dimension(3,(im+1)*(jm+1)) :: xi,xj
      Integer i,j,n
      Character*100 fname,fnext


      write(unit=fnext,fmt=100) 1000+irank
 100  Format(i4)
      fname = 'TEST_'//trim(fnext)//'.dat'
      if( irank.eq.0 ) PRINT*,'filename:',trim(fname)

      open(unit=10,file=fname,status='unknown',form='formatted')
      write(10,*) 'variables = x y z'
      write(10,*) 'zone T="',irank,'",',
     &            'N=',(im+1)*(jm+1),',',
     &            'E=',(im  )*(jm  ),',',
     &            'F=FEPOINT, ET=QUADRILATERAL'

      Do j = 1,jm+1
      Do i = 1,im+1
         n = (j-1)*(im+1) + i
         write(10,*) xi(:,n)
      Enddo
      Enddo
      Do n = 1,im*jm
         write(10,*) icon(1:4,n)
      Enddo

      write(10,*) 'zone T="',irank,'",',
     &            'N=',(im+1)*(jm+1),',',
     &            'E=',(im  )*(jm  )*2,',',
     &            'F=FEPOINT, ET=TRIANGLE'

      Do j = 1,jm+1
      Do i = 1,im+1
         n = (j-1)*(im+1) + i
         write(10,*) xj(:,n)
      Enddo
      Enddo

      Do n = 1,im*jm*2
         write(10,*) jcon(1:3,n)
      Enddo

      close(10)

      return
      end

      Subroutine match_test(
     &                 icomw, irank, nproc, im,jm, icon, xi, jcon, xj,
     &                 num_recv, num_send, irecv, isend, area,
     &                 irdis, ircnt, isdis, iscnt )
      use MPI
      Implicit none
      Integer icomw, irank, nproc, im,jm
      Integer,dimension(0:4,im*jm) :: icon
      Integer,dimension(0:4,im*jm*2) :: jcon
      Real*8, dimension(3,(im+1)*(jm+1)) :: xi,xj
      Integer num_recv, num_send
      Integer,dimension(0:nproc-1) :: irdis,ircnt, isdis,iscnt
      Integer,dimension(num_recv) :: irecv
      Real*8, dimension(num_recv) :: area
      Integer,dimension(num_send) :: isend
      Real*8, allocatable, dimension(:) :: rdum,sdum,q,p
      Integer,allocatable, dimension(:) :: ireq
      Integer i,n, ier
      Character*100 fname,fnext


      allocate( q(im*jm*2) )
      allocate( p(im*jm) )
      allocate( rdum(num_recv) )
      allocate( sdum(num_send) )
      allocate( ireq(0:2*nproc-1) )

      do n = 1,im*jm*2
      !  q(n) = xj(1, jcon(1,n) )    ! test to be completed
         q(n) = dble(irank)  ! HACK
      enddo

      do n = 0,nproc-1
         do i = isdis(n)+1,isdis(n)+iscnt(n)
            sdum(i) = q( isend(i)+1 )
         enddo
      enddo

      rdum(:) = -9999999.0d0    ! a dummy number to trap errors
      do n = 0,nproc-1
         call MPI_Irecv( rdum( irdis(n)+1 ), ircnt(n),
     &                   MPI_DOUBLE_PRECISION,n,n+1000,
     &                   icomw, ireq(n), ier )
      enddo
      do n = 0,nproc-1
         call MPI_Issend( sdum( isdis(n)+1 ), iscnt(n),
     &                   MPI_DOUBLE_PRECISION,n,irank+1000,
     &                    icomw, ireq(n+nproc), ier )
      enddo
      n = 2*nproc
      call MPI_WAITALL( n, ireq, MPI_STATUSES_IGNORE, ier )

      deallocate( ireq )

      p(:) = 0.0d0    ! initialization of accumulating sum
      do n = 0,nproc-1
         do i = irdis(n)+1,irdis(n)+ircnt(n)
            p( irecv(i)+1 ) = p( irecv(i)+1 ) + area(i)*rdum(i)
         enddo
      enddo

      write(unit=fnext,fmt=100) 1000+irank
 100  Format(i4)
      fname = 'POST_'//trim(fnext)//'.dat'
      if( irank.eq.0 ) PRINT*,'filename:',trim(fname)

      open(unit=10,file=fname,status='unknown',form='formatted')
      write(10,*) 'variables = x y z val'
      write(10,*) 'zone T="got_',irank,'",',
     &            'N=',(im+1)*(jm+1),',',
     &            'E=',(im  )*(jm  ),',',
     &            'VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)',
     &            'ZONETYPE=FEQUADRILATERAL, DATAPACKING=BLOCK'

      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(1,n)
      Enddo
      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(2,n)
      Enddo
      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(3,n)
      Enddo
      Do n = 1,(jm  )*(im  )
         write(10,*) p(n)
      Enddo

      Do n = 1,im*jm
         write(10,*) icon(1:4,n)
      Enddo

      write(10,*) 'zone T="impose_',irank,'",',
     &            'N=',(im+1)*(jm+1),',',
     &            'E=',(im  )*(jm  )*2,',',
     &            'VARLOCATION=([1-3]=NODAL,[4]=CELLCENTERED)',
     &            'ZONETYPE=FETRIANGLE, DATAPACKING=BLOCK'

      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(1,n)
      Enddo
      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(2,n)
      Enddo
      Do n = 1,(jm+1)*(im+1)
         write(10,*) xi(3,n)
      Enddo
      Do n = 1,(jm  )*(im  )*2
         write(10,*) q(n)
      Enddo

      Do n = 1,im*jm*2
         write(10,*) jcon(1:3,n)
      Enddo

      deallocate( q )
      deallocate( p )
      deallocate( sdum )
      deallocate( rdum )

      close(10)

      return
      end

