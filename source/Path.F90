!
!
!****************************************************************************
!
!  SUBROUTINE: Path
!
!  PURPOSE:  Compute path of the flux.
!
!****************************************************************************
SUBROUTINE Path(ixini,iyini,acumulado,acumuladoVel,velocidades,control,topostat,topo)
!
!
!	Libreria FORTRAN, necesaria para el generador de numeros pseudoaleatorios
	!USE DFLIB
!
!	Variables globales
	use DebrisDiceGlobals

!
!	Las variables reales son de doble precisión e implicitas
    implicit double precision (a-h,o-z)
!
!	Variables explicitas
	real rand			!Numero REAL*4
!
!	Variables dinamicas (tamaño variable)
    dimension topo(mx,my)					!Malla de topografia
    dimension topostat(mx,my,8)				!Malla de pesos
	dimension acumulado(mx,my)				!Malla de resultados
	dimension acumuladoVel(mx,my)			!Malla de resultados velocidad
	dimension velocidades(mx,my)			!Malla de velocidades
	LOGICAL, dimension(mx,my) :: control	!Malla de control de borde y sumideros (BOOLEAN)
	LOGICAL :: controlVel					!Calculo de velocidades
	LOGICAL :: controlAng					!Calculo de velocidades
!
!
!	Bucle principal de cálculo
	do i = 1, iter
!
!
!		Posiciones iniciales
		ix = ixini
		iy = iyini
!
!		Variables utilizadas para angle abast
		dist = 0.d0
		cota = 0.d0
		controlAng = .TRUE.
!
!		Variables para velocidad
		velCuad = 0.d0
		controlVel = .TRUE.
!
!		Actualizamos los resultados para incluir el punto de inicio
		acumulado(ix,iy) = acumulado(ix,iy) + 1.d0
		acumuladoVel(ix,iy) = acumuladoVel(ix,iy) + 1.d0
!
!		Generar un número aleatorio no es facil, para ello este FORTRAN implementa el algoritmo:
!			Prime Modulus M Multiplicative Linear Congruential Generator
!			Park and Miller, "Random Number Generators: Good Ones Are Hard to Find," CACM, October 1988, Vol. 31, No. 10.) 
!
!		Primero generamos la semilla, constante para cada iteración Montecarlo		
		CALL RANDOM_SEED()
!
!		Continuamos hasta salir
		LOOP_PATH:do while (controlAng .OR. controlVel)
!
!			Cota actual
			zini = topo(ix,iy)
!
!			Numero aleatorio (rand)		
			CALL RANDOM_NUMBER(rand)
			ran = DBLE(rand)	!Convertimos a REAL*8
!
!			Determinamos la dirección del flujo
			LOOP_SEARCH:do k= 1,8
!			
				if (topostat(ix,iy,k) .GE. ran) EXIT LOOP_SEARCH
!	
!			Final loop search
			enddo LOOP_SEARCH
!
!			Determinamos la celda de destino del flujo
			CHECK_PATH: SELECT CASE (k)

				CASE (1)
					iy = iy + 1
					auxdist = dy

				CASE (2)
					ix = ix + 1
					iy = iy + 1
					auxdist = dd
				
				CASE (3)
					ix = ix + 1
					auxdist = dx

				CASE (4)
					ix = ix + 1
					iy = iy - 1
					auxdist = dd

				CASE (5)
					iy = iy - 1
					auxdist = dy

				CASE (6)
					ix = ix - 1
					iy = iy - 1
					auxdist = dd
				
				CASE (7)
					ix = ix - 1
					auxdist = dx

				CASE (8)
					ix = ix - 1
					iy = iy + 1
					auxdist = dd

			END SELECT CHECK_PATH
!
!			Actualizamos longitud
			dist = dist + auxdist
!
!			Actualizamos los el estado
			if (.NOT. control(ix,iy)) then
				controlAng = .FALSE.
				controlVel = .FALSE.
			endif
!
!			Si seguimos dentro sumamos
			if (controlAng) then
				acumulado(ix,iy) = acumulado(ix,iy) + 1.d0
			endif
!
!			Controlamos angle abast
			cota = topo(ixini,iyini) - topo(ix,iy)
			angleAbast = cota / dist
!
!			La caida debe ser mayor que l'abast
!			if (angleAbast .LT. abast) EXIT LOOP_PATH
			if (angleAbast .LT. abast) then
				controlAng = .FALSE.
			endif
!
!			Calculo de velocidades
			if (controlVel) then
!
!				Cota final
				zfin = topo(ix,iy)
				tanTheta = (zini - zfin) / auxdist
				sintheta = DSIN(DATAN(tanTheta))
				cosTheta = DCOS(DATAN(tanTheta))
				if (iFricType .EQ. 1) then
					velCuad = velCuad + 2.d0 * auxdist * (9.81d0 * (sintheta - rmu*cosTheta) - DABS(velCuad) / rk)
				elseif(iFricType .EQ. 2) then
					Friccion =  shBG * (1.d0 + cte_BG * (DSQRT(DABS(velCuad)))**expm)
					velCuad = velCuad + auxdist * (9.81d0 * sintheta - Friccion / bingDepth)
				else
					velCuad = velCuad + auxdist * (9.81d0 * sintheta - viscBG * 3.d0 / bingDepth**2)
				endif
!
!				Solo velocidades positivas
				if (velCuad .GT. 0.d0) then
					velocidades(ix,iy) = velocidades(ix,iy) + velCuad
					acumuladoVel(ix,iy) = acumuladoVel(ix,iy) + 1.d0
				else
					controlVel = .FALSE.
				endif
			endif
!
!
		enddo LOOP_PATH	!Final DO WHILE
!
	enddo	!Final DO iter
!
!
end
!
!
!