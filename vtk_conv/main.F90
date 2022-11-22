program main
  use omp_lib 
  implicit none

  integer :: n,i,nthreads,np,npp,ngrab,Nframes,iframe,i_loop_finish,N_start,ii,i_PART_counter
  integer,parameter :: np_max = 999999
  integer, parameter :: i_PART_counter_max=20000
  character chartemp*40, name_orig*40
  character name_vtu*40, name_vtu2*12, name_vtu3*9
  character np_string3*3, np_string4*4, np_string5*5
  character np_string6*6, np_string7*7, np_string8*8
  character frame_string1*1, frame_string2*2, frame_string3*3
  character frame_string4*4, frame_string5*5, frame_string6*6
  character supp*4,supp3*3,supp2*2,supp1*1, zero_string
  character string1*100,string2*100,string3*100,string4*100
  character chartemp2*100
  character proc5*5
  CHARACTER(LEN=10) :: FMT,FMT1
  CHARACTER(LEN=1)  :: TAB,DQ
      
  integer :: itn,ifi,ifo,di
  real :: dr
      
  real,allocatable,dimension(:):: xp,up,vp,vort,Temp,yp,alpha
  real time(i_PART_counter_max), DT(i_PART_counter_max)
  integer np_all(i_PART_counter_max), IT(i_PART_counter_max)
  integer processor(np_max),node_type(np_max),subset_flag
  real  DT1(i_PART_counter_max),DT2(i_PART_counter_max)  
  integer :: nprocs,iproc,np_ini,np_end
  
  allocate(xp(np_max))
  allocate(yp(np_max))
  allocate(up(np_max))
  allocate(vp(np_max))
  yp = 0.0d0
  allocate(vort(np_max))
  allocate(Temp(np_max))
  allocate(alpha(np_max))

  TAB=CHAR(9)     
  FMT="(A)"
  FMT1="(2A)"
  DQ=CHAR(34)
  

 
      
     
  !! LOAD AND READ TIME,DT FROM FILE DT. 
  open(unit=70,file='../data_out/time_out',status='old')
  i_loop_finish = 0
  i_PART_counter = 0
  do while(i_loop_finish.eq.0)
     i_PART_counter = i_PART_counter + 1
     if(i_PART_counter.gt.i_PART_counter_max)then
        write(6,*) 'Number of entries in file DT exceeds max value'
        write(6,*) 'i_PART_counter.gt.i_PART_counter_max'
        write(6,*) 'Adjust i_PART_counter_max, i_PART_counter_max = ',i_PART_counter_max
        stop
     endif
     read(70,*,END = 76)time(i_PART_counter),np_all(i_PART_counter),IT(i_PART_counter),DT(i_PART_counter),nprocs
           
     !Determine whether to exit loop
     if(i_loop_finish.eq.0)then
        i_loop_finish = i_loop_finish - 1
     endif
76   i_loop_finish = i_loop_finish + 1
  end do
  N_start = 1
  Nframes = i_PART_counter-2  !Why -2?
  write(6,*) "There are ",Nframes+1,"frames."

  !! JRCK addition...          
  write(6,*) "Enter starting frame"
  read(*,*) N_start
    
  ngrab = N_start-1
        
  !! Loop over each frame   
!  !$omp parallel do private(name_vtu,ngrab,npp,iproc,name_orig,supp1,supp2,supp3,supp, &
!  !$omp xp,zp,h,node_type,ro,up,vp,vort,energy,Temp,Y0,np_ini,np_end,np,i,processor, &
!  !$omp string1,string2,string3,string4,np_string3,np_string4,np_string5,np_string6,np_string7,np_string8, &
!  !$omp itn,ifi,ifo,proc5)
  do iframe=N_start,Nframes+1

     !! Each thread needs a different io number for in and out files
     itn = 0!omp_get_thread_num()
     ifi = 23+itn
     ifo = 123+itn

     write(supp,'(i4.4)') ngrab
     name_vtu ='../paraview_files/LAYER'//supp//'.vtu'
     open(ifo,file=name_vtu,status='unknown')

     !! Increment the frame counter ngrab
     ngrab=iframe
     npp=0 ! keeps track of total number of particles/nodes across all processors
     !! Loop over each processor
             
!       % READ IN THE PART FILE FOR EACH FRAME
        if(ngrab.lt.10) then
           write(supp1,'(i0)') ngrab
           name_orig='../data_out/layer_'//supp1
        end if
        if(ngrab.ge.10.and.ngrab.lt.100) then
           write(supp2,'(i2)') ngrab
           name_orig='../data_out/layer_'//supp2
        end if
        if(ngrab.ge.100.and.ngrab.lt.1000) then
           write(supp3,'(i3)') ngrab
           name_orig='../data_out/layer_'//supp3
        end if
        if(ngrab.ge.1000.and.ngrab.lt.10000) then
           write(supp,'(i4)') ngrab
           name_orig='../data_out/layer_'//supp
        end if
        
        open(ifi,file=name_orig,status='old')
!       % READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES                      
        read(ifi,*) np
        np_ini = npp + 1
        np_end = np_ini + np 
           do i=np_ini,np_end
              read(ifi,*,end=300) xp(i),yp(i),up(i),vp(i),vort(i),Temp(i),alpha(i)
              npp=npp+1
           enddo
300     close (ifi)


     np = npp

     write(6,*) "Frame",iframe,"with ",np,"particles, from ",nprocs,"processors."
                                                                      
201  format(a40)
202  format(a100)
203  format(a25,i7,a17,i7,a2)
211  format(a21)
!     % OUTPUT TO FILE IN VTU FORMAT 
     if(np.lt.1000)then       
        write(np_string3,'(i3.3)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string3//DQ//' NumberOfCells='//DQ//np_string3//DQ//'>'
     elseif(np.lt.10000)then       
        write(np_string4,'(i4.4)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string4//DQ//' NumberOfCells='//DQ//np_string4//DQ//'>'
     elseif(np.lt.100000)then       
        write(np_string5,'(i5.5)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string5//DQ//' NumberOfCells='//DQ//np_string5//DQ//'>'
     elseif(np.lt.1000000)then       
        write(np_string6,'(i6.6)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string6//DQ//' NumberOfCells='//DQ//np_string6//DQ//'>'
     elseif(np.lt.10000000)then       
        write(np_string7,'(i7.7)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string7//DQ//' NumberOfCells='//DQ//np_string7//DQ//'>'
     elseif(np.lt.100000000)then       
        write(np_string8,'(i8.8)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string8//DQ//' NumberOfCells='//DQ//np_string8//DQ//'>'
     else
        write(6,*) 'Too many particles for np_string'
        stop  
     endif
 
     string1 = '<?xml version='//DQ//'1.0'//DQ//'?>'
     string2 = '<VTKFile type= '//DQ//'UnstructuredGrid'//DQ//'  version= '//DQ//'0.1'//DQ//&
               '  byte_order= '//DQ//'BigEndian'//DQ//'>'
     string3 = ' <UnstructuredGrid>'
     write(ifo,211)string1
     write(ifo,202)string2
     write(ifo,202)string3
     write(ifo,202)string4
              
     !! Start of point data
     string1 = '   <PointData Scalars='//DQ//'Pressure'//DQ//' Vectors='//DQ//'Velocity'//DQ//'>'
     write(ifo,202)string1

     string2 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'h'//DQ//' format='//DQ//'ascii'//DQ//'>'

     !! Vorticity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Vorticity'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)vort(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3


     !! Temperature
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Temperature'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)Temp(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
     
     !! spare output
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'alpha'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)alpha(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3     

     !! U-velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'u'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)up(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3

     !! V-velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'v'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202)string1
     do ii=1,np
        write(ifo,*)vp(ii)
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
             
     !! Vector velocity
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//DQ//'Velocity'//DQ// &
               ' NumberOfComponents='//DQ//'3'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)up(ii),vp(ii),0.0d0
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3


     !! End of point data
     string4 = '   </PointData>'
     write(ifo,202) string4
              
     !! Finally, particle positions!!
     string2 = '   <Points>'
     string1 = '    <DataArray type='//DQ//'Float32'//DQ//' NumberOfComponents='//DQ//'3'//DQ// &
               ' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string2
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)xp(ii),yp(ii),0.0d0
     enddo
     string3 = '    </DataArray>'
     string2 = '   </Points>'
     write(ifo,202) string3
     write(ifo,202) string2

     !! WRITE CELL DATA. CELL IS OF TYPE VERTEX.          
     string2 = '   <Cells>'
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'connectivity'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string2
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)ii-1
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
       
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'offsets'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)ii
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
             
     string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ//'types'//DQ//' format='//DQ//'ascii'//DQ//'>'
     write(ifo,202) string1
     do ii=1,np
        write(ifo,*)1
     enddo
     string3 = '    </DataArray>'
     write(ifo,202) string3
        
     !! Final bits      
     string1 = '   </Cells>' 
     string2 = '  </Piece>'
     string3 = ' </UnstructuredGrid>'
     string4 = '</VTKFile>'
     write(ifo,202) string1
     write(ifo,202) string2
     write(ifo,202) string3
     write(ifo,202) string4
     close(24)


  enddo
!  !$omp end parallel do



  deallocate(xp,yp,up,vp)
  deallocate(vort,Temp)  


  stop
end program main
!! ------------------------------------------------------------------------------------------------

