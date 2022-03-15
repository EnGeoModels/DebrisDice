!
!
!****************************************************************************
!
!  SUBROUTINE: StatTopo
!
!  PURPOSE:  Compute topostat.
!
!****************************************************************************
SUBROUTINE StatTopo(topostat,topo)
!
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
    dimension topostat(mx,my,8)				!Malla de pesos
!
!
!
!	Determinamos los 8 ángulos de contorno
    do i = 2, (mx - 1)
		do j = 2, (my - 1)
			do k = 1, 8
!
!			El ángulo lo podemos determinar de 2 maneras, con la TAN o con el SEN
	CHECK_ANGLE: SELECT CASE (k)

				CASE (1)
!
!					Seno
					!delta_H = topo(i,j+1)   - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dy**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente					
					topostat(i,j,k) = (topo(i,j+1)   - topo(i,j)) / dy

				CASE (2)
!
!					Seno
					!delta_H = topo(i+1,j+1) - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dd**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i+1,j+1) - topo(i,j)) / dd
				
				CASE (3)
!
!					Seno
					!delta_H = topo(i+1,j)   - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dx**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i+1,j)   - topo(i,j)) / dx

				CASE (4)
!
!					Seno
					!delta_H = topo(i+1,j-1) - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dd**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i+1,j-1) - topo(i,j)) / dd

				CASE (5)
!
!					Seno
					!delta_H = topo(i,j-1)   - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dy**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i,j-1)   - topo(i,j)) / dy

				CASE (6)
!
!					Seno
					!delta_H = topo(i-1,j-1) - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dd**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i-1,j-1) - topo(i,j)) / dd
				
				CASE (7)
!
!					Seno
					!delta_H = topo(i-1,j)   - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dx**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i-1,j)   - topo(i,j)) / dx

				CASE (8)
!
!					Seno
					!delta_H = topo(i-1,j+1) - topo(i,j)
					!HIP = DSQRT(delta_H**2 + dd**2)
					!topostat(i,j,k) = delta_H / HIP
!
!					Tangente
					topostat(i,j,k) = (topo(i-1,j+1) - topo(i,j)) / dd

    END SELECT CHECK_ANGLE
!
!				Corrección para fondo plano
				if (topostat(i,j,k) .EQ. 0.d0) topostat(i,j,k) = -0.001
!
!				Corregimos para evitar upward
				topostat(i,j,k) = -DMIN1(topostat(i,j,k),0.d0)
!
!				Ampliamos la dispersion
				!topostat(i,j,k) = DSQRT(topostat(i,j,k))


			enddo
		enddo
    end do
!
!
end
!
!
!
!
!
