!
!****************************************************************************
!
!  SUBROUTINE: Normalize
!
!  PURPOSE:  Normalize topostat.
!
!****************************************************************************
SUBROUTINE Normalize(topostat,control,topo)
!
!
!	Variables globales
	use DebrisDiceGlobals
!
!
	implicit double precision (a-h,o-z)
!
!
!	Variables dinamicas (tamaño variable)
    dimension topostat(mx,my,8)				!Malla de pesos
    dimension topo(mx,my)					!Malla de topografia
	LOGICAL, dimension(mx,my) :: control	!Malla de control de borde y sumideros (BOOLEAN)
!
!
!
!	Determinamos los coeficientes
!	Para ello utilizamos la probabilidad acumulada
    do i = 2, (mx - 1)
		do j = 2, (my - 1)
!
!			Calculamos los pesos normalizando
			suma = SUM(topostat(i,j,:))
			amax = MAXVAL(topostat(i,j,:))
!
			if ( (amax .GT. amin) .AND. (topo(i,j) .NE. DBLE(nodata)) ) then
				topostat(i,j,:) = topostat(i,j,:) / suma
				control(i,j)    = .TRUE.
			else
				topostat(i,j,:) = 0.d0
				control(i,j)	= .FALSE.
			endif
!
!			Calculamos probabilidad acumulada
			do k = 2, 8
				topostat(i,j,k) = topostat(i,j,k) + topostat(i,j,k-1)
			enddo
!
		enddo
	enddo
!
!	El borde de la malla es control de salida
	control(:,1)   = .FALSE.
	control(:,my)  = .FALSE.
	control(1,:)   = .FALSE.
	control(mx,:)  = .FALSE.
!
!
!
end
!
