!Subrutina que escribe los resultados en forma de GRID Arcview
!para poder visualizarlos

!==================================================
subroutine GridOut(fname, ArrayToWrite, mx, my, dx, xcorner, ycorner)
!==================================================
!
!
	implicit double precision (a-h,o-z)
!
!	Local variables
	REAL*8,   DIMENSION(mx,my) :: ArrayToWrite                !Matriz donde se almacenan las variables
	character*30 fname
!
		open(unit=100,file=fname,status='unknown',form='formatted')
!
!		Keywords de la GRID
		write(100,1000) 'ncols         ', mx
		write(100,1000) 'nrows         ', my  
!
		write(100,1001) 'xllcorner     ', xcorner
		write(100,1001) 'yllcorner     ', ycorner
!    
		write(100,1001) 'cellsize      ', dx
		write(100,1000) 'NODATA_value  ', -9999
!
!		Escribimos la malla
		do j = 1,my
			write(100,1002) (ArrayToWrite(i,j), i =1,mx-1)
			write(100,1004) ArrayToWrite(mx,j)
		end do
!    
!		Cerramos fichero
		close(100)

!
!	Formatos de Keyword
1000 format(A14, I10)
1001 format(A14, F14.6)
1002 format(F12.7, $)
1003 format(/)    
1004 format(F12.7)
!
!
end subroutine

