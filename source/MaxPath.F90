!
!
!****************************************************************************
!
!  SUBROUTINE: MaxPath
!
!  PURPOSE:  Compute max slope path of the flux.
!
!****************************************************************************
SUBROUTINE MaxPath(ixini,iyini,acumulado,control,topo)
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
    dimension topo(mx,my)					!Malla de topografia
	dimension acumulado(mx,my)				!Malla de resultados
	LOGICAL, dimension(mx,my) :: control	!Malla de control de borde y sumideros (BOOLEAN)
!
!	Posiciones iniciales
	i = ixini
	j = iyini
!
!	Actualizamos los resultados para incluir el punto de inicio
	acumulado(i,j) = 1.d0
!
!	Continuamos hasta salir
	LOOP_PATH:do while (control(i,j))
!
!		Direccion de máxima pendiente
		dir = 0.d0
		inode = 0
!
!		Determinamos la dirección del flujo
		LOOP_SEARCH:do k= 1,8
!			
!			Determinamos la máxima pendiente
			CHECK_ANGLE: SELECT CASE (k)

				CASE (1)
!
!					Tangente					
					topostat = (topo(i,j+1)   - topo(i,j)) / dy

				CASE (2)
!
!					Tangente
					topostat = (topo(i+1,j+1) - topo(i,j)) / dd
				
				CASE (3)
!
!					Tangente
					topostat = (topo(i+1,j)   - topo(i,j)) / dx

				CASE (4)
!
!					Tangente
					topostat = (topo(i+1,j-1) - topo(i,j)) / dd

				CASE (5)
!
!					Tangente
					topostat = (topo(i,j-1)   - topo(i,j)) / dy

				CASE (6)
!
!					Tangente
					topostat = (topo(i-1,j-1) - topo(i,j)) / dd
				
				CASE (7)
!
!					Tangente
					topostat = (topo(i-1,j)   - topo(i,j)) / dx

				CASE (8)
!
!					Tangente
					topostat = (topo(i-1,j+1) - topo(i,j)) / dd

			END SELECT CHECK_ANGLE
!
!			Selecionamos el máximo
			if(topostat .LT. dir) then
				dir = topostat
				inode = k
			endif
!	
!		Final loop search
		enddo LOOP_SEARCH
!
!		Solo bajamos
		if (inode .EQ. 0) EXIT LOOP_PATH
!
!		Determinamos la celda de destino del flujo
		CHECK_PATH: SELECT CASE (inode)

			CASE (1)
				j = j + 1

			CASE (2)
				i = i + 1
				j = j + 1
				
			CASE (3)
				i = i + 1

			CASE (4)
				i = i + 1
				j = j - 1

			CASE (5)
				j = j - 1

			CASE (6)
				i = i - 1
				j = j - 1
				
			CASE (7)
				i = i - 1

			CASE (8)
				i = i - 1
				j = j + 1

		END SELECT CHECK_PATH
!
!		Actualizamos los resultados
		acumulado(i,j) = 1.d0
!
!
	enddo LOOP_PATH	!Final DO WHILE
!
!
end
!
!
