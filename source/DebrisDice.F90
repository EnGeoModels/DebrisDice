!  DebrisDice.f90 
!
!  FUNCTIONS:
!	DebrisDice      - Entry point of console application.
!
!
!	Leemos la topografia
!	Leemos el fichero de datos
!	Creamos los pesos estadísticos
!	Calculamos los casos
!
!	La matriz acumulado(:,:) es el resultado (tipo DOUBLE)
!	La matriz control(:,:) detecta el borde y sumideros (pozos) (tipo BOOLEAN)
!	

!****************************************************************************
!
!  PROGRAM: DebrisDice
!
!  PURPOSE:  Compute debris flow using statistical aproach.
!
!****************************************************************************

program DebrisDice
!
!	Variables globales
	use DebrisDiceGlobals
!
	implicit double precision (a-h,o-z)
!
!
!	Variables dinamicas (tamaño variable)
    ALLOCATABLE :: topo(:,:)					!Malla de topografia
    ALLOCATABLE :: topostat(:,:,:)				!Malla de pesos
	ALLOCATABLE :: acumulado(:,:)				!Malla de resultados
	ALLOCATABLE :: acumuladoVel(:,:)			!Malla de resultados con velocidad
	ALLOCATABLE :: isources(:,:)				!Malla de puntos de inicio
	ALLOCATABLE :: velocidades(:,:)				!Malla de velocidades
	LOGICAL,   ALLOCATABLE :: control(:,:)		!Malla de control de borde y sumideros (BOOLEAN)
!
!	Cabecera
	write(6,'("**************************************************")')
	write(6,'("        DebrisDice method for hazard analisys     ")')
	write(6,'("**************************************************",/)')
!
!
!
!	Fichero de topografia
    open(100,file='topo.dat',status='old',ACTION='read',IOSTAT=istat)
!
!   Exception handling
    if (istat .ne. 0) then
        write(*,*) "topo.dat cannot be opened"
        stop 1
    end if    
!
	write(6,'("Reading topo data and constructing matix",/)')
!
!	Entrada de la topografia topo.dat
    read(100,*) dummy,mx        !numero nodos X
    read(100,*) dummy,my        !numero nodos Y
    read(100,*) dummy,xcorner   !Coordenada X esquina superior izq.
    read(100,*) dummy,ycorner   !Coordenada X esquina superior izq.
    read(100,*) dummy,dx        !Delta X
    read(100,*) dummy,nodata	!nodata ESRI
!
!	Calculamos el resto de distandias
	dy = dx
	dd = DSQRT(2.d0) * dx	!Distancia diagonal
!
!
!	Dimensionamos las variables dinámicas
	ALLOCATE(topo(mx,my))			!Malla de topografia
	ALLOCATE(topostat(mx,my,8))		!Malla de pesos
	ALLOCATE(acumulado(mx,my))		!Malla de resultados
	ALLOCATE(acumuladoVel(mx,my))	!Malla de resultados con velocidad
	ALLOCATE(control(mx,my))		!Malla de control de borde y sumideros
	ALLOCATE(isources(mx,my))		!Malla de puntos de inicio
	ALLOCATE(velocidades(mx,my))	!Malla de velocidades
!
!	Leemos la malla de terreno
    do j = 1, my
        read(100,*) (topo(i,j), i =1, mx)
    end do
!   
!	Cerramos el fichero
    close(100)   
!
!
!
!	Abrimos ficheros de entrada de datos de control
	open(100,file='input.dat',status='old',form='formatted',ACTION='read',IOSTAT=istat)
!
!   Exception handling
    if (istat .ne. 0) then
        write(*,*) "input.dat cannot be opened"
        stop 1
    end if
!
	write(6,'("Reading configuration file...",/)')
!
!	Entrada de parametros de input.dat:
	read(100,*) itipo,imaxpend		!Calculo para un punto o para una malla, calcular máxima pendiente
	write(6,'("Type:               ",I5," Steepest path: ",I5)') itipo,imaxpend
    read(100,*) iter				!Numero de iteraciones
	write(6,'("Iterations:         ",I5)') iter
!
!	Selecionamos si se trata de un calculo para un punto o una malla
	if (itipo .EQ. 1) then
		read(100,*) rxini,ryini		!Coordenadas UTM x e y de la celda de origen 
		write(6,'("X:       ",F12.3," Y: ",F12.3)') rxini,ryini
    else
		read(100,*) fname			!Nombre del fichero de la malla de puntos de inicio	
		write(6,'("File: ",A20)') fname
	endif
!
	read(100,*) abast, amin			!Maxim abast, Angulo minimo de parada
	write(6,'("Range angle: ",F12.3," Min slope: ",F12.3)') abast, amin
!
!	Tipo friccion
	read(100,*) iFricType			!Frictional or viscous
	write(6,'("Friction type:      ",I5)') iFricType
!
!	Voellmy
	read(100,*) rmu, rk				!Coulomb friction factor, mass to drag ratio
	write(6,'("Mu:          ",F12.3," Drag to mass: ",F12.3)') rmu, rk
!
!	Herschel-Bulkley
	read(100,*) shBG		!ShearYield, damos la tension en Pa
	write(6,'("ShearYield:  ",F12.3)') shBG

	read(100,*) viscBG		!Viscosidad, se da como dinamica (Pa·s)
	write(6,'("Viscosity:   ",F12.3)') viscBG

	read(100,*) expBG		!Exponente Herschel-Bulkley
	write(6,'("Exponent:    ",F12.3)') expBG

	read(100,*) densS		!Density (kg/m3)
	write(6,'("Density:     ",F12.3)') densS

	read(100,*) bingDepth	!Depth used in viscous velocity calculation
	write(6,'("Depth:       ",F12.3)') bingDepth

!
!	Cerramos el fichero
    close(100)
!
!	Dividimos por la densidad
	viscBG = viscBG / densS
	shBG   = shBG   / densS
!
!
!	Seleccionamos en funcion de si calculamos un source o una malla de sources
	if (itipo .EQ. 1) then
!
!		Comprobamos que las coordenadas UTM esten dentro de la malla
		rmaxUTMx = xcorner + DBLE(mx) * dx
		rmaxUTMy = ycorner + DBLE(my) * dy
!
		if (rxini .GE. rmaxUTMx .OR. ryini .GE. rmaxUTMy) then
			write(6,'("ERROR: UTM out of grid domain",/)')
			goto 1000
		endif
!
!		Convertimos a coordenadas UTM a numero de celda
		ixini = IDINT((rxini - xcorner) / dx)
		iyini = IDINT((ryini - ycorner) / dy)
!
!		Las coordenadas del nodo origen estan en referencia dextrogira
!		pasamos a levogira (las mallas ARCVIEW tienen y=1 arriba e y=my abajo, hay que invertir
		iyini = my - iyini + 1
	else
!
!		Malla de sources
		open(100,file=fname,status='old')
!
		write(6,'("Reading init points grid...",/)')
!
!		Comprobamos que sea el mismo header que la topografia
		read(100,*) dummy,imx			!numero nodos X
		if(imx.ne.mx) goto 1000
    
		read(100,*) dummy,imy			!numero nodos Y
		if(imy.ne.my) goto 1000
    
		read(100,*) dummy,txcorner		!Coordenada X esquina superior izq.
		if(DABS(txcorner - xcorner) .GT. 1.d0) goto 1000
    
		read(100,*) dummy,tycorner      !Coordenada Y esquina superior izq.
		if(DABS(tycorner - ycorner) .GT. 1.d0) goto 1000
    
		read(100,*) dummy,tdx           !Delta X
		if(DABS(tdx - dx) .GT. 1.d-3) goto 1000

		read(100,*) dummy,nodata        !nodata ESRI
!
!		Leemos la malla de terreno
		do j = 1, my
			read(100,*) (isources(i,j), i =1, mx)
		end do
!   
!		Cerramos el fichero
		close(100)   
!
	endif
!
!
!
!	Preprocesamos la geometría
!
	write(6,'("Geometry preprocessing...",/)')
!
	Call StatTopo(topostat,topo)
!
	write(6,'("Statistics computation...",/)')
!
	Call Normalize(topostat,control,topo)
!
!
!
!	Inicializamos los resultados
	acumulado(:,:) = 0.d0
	acumuladoVel(:,:) = 0.d0
	velocidades(:,:) = 0.d0
!
!	Valores reologicos relacionados con BINGHAM, son independientes del nodo
	if(iFricType .EQ. 2) then
!
!		Tensión umbral no nula
		if(shBG .GT. 0.d0) then
			rm   = 1.d0 / expBG
			expm = 1.d0 / (rm + 0.15d0)
			ctem = ((rm + 1.d0)*(rm + 2.d0) / ((0.74d0 + 0.656d0 * rm)*bingDepth)) ** expm
			cocm = (viscBG / shBG)**(expm*rm)
			cte_BG = ctem * cocm
!
!		Tension umbral nula --> Newtoniano!
		else
			iFricType = 3
		endif
	endif
!
!     
	write(6,*) 'Iterations init'
!
!
!	Diferenciamos el tipo de calculo
	if (itipo .EQ. 1) then
!
!		Un solo source
		write(6,*) 'One source.'
		Call Path(ixini,iyini,acumulado,acumuladoVel,velocidades,control,topostat,topo)	
!
    else
!
!		Multiples sources
		totalNumber = SUM(MAX(isources,0))
		iNumber = 0
!
		do i = 1,mx
			do j = 1,my
!
!				Si se trata de u source calaculamos
				if(isources(i,j) .EQ. 1) then
					iNumber = iNumber + 1
					write(6,'("Done:     ",F12.3)') (DBLE(iNumber) / DBLE(totalNumber) * 100.d0)			
!
					ixini = i
					iyini = j
					Call Path(ixini,iyini,acumulado,acumuladoVel,velocidades,control,topostat,topo)	
				endif
!
			enddo
		enddo

	endif
!
!
!
!	Salida de resultados mediante GRID Arcview
	write(6,'("Topo output",/)')
!

!	Escribimos resultados
    fname = 'Topo.txt'
!
!	Escribimos la malla
	call GridOut(fname, topo, mx, my, dx, xcorner, ycorner)
!
!
!	Salida de resultados velocidad mediante GRID Arcview
	write(6,'("Velocity result output",/)')

!   
!	Escribimos resultados
    fname = 'Velocity.txt'
!
!	Calculamos la velocidad media
	do j = 1, my
		do i = 1, mx
!
			velCuad = velocidades(i,j)
			acum = DMAX1(acumuladoVel(i,j),1.d0)
!
			if(velCuad .GT. 0.d0) then
				velocidades(i,j) = DSQRT(velCuad / acum)
			else
				velocidades(i,j) = -DSQRT(-velCuad / acum)
			endif
!
		enddo
	end do
!
!	Escribimos la malla
	call GridOut(fname, velocidades, mx, my, dx, xcorner, ycorner)
!
!
!	Salida de resultados mediante GRID Arcview
	write(6,'("Result output",/)')

!   
!	Escribimos resultados
    fname = 'Result_Angle.txt'
!
!	Normalizamos los resultados
	rmaxval = DMAX1(MAXVAL(acumulado),1.d0)
	acumulado = acumulado / rmaxval
!
!	Escribimos la malla
	call GridOut(fname, acumulado, mx, my, dx, xcorner, ycorner)
!
!
!	Salida de resultados mediante GRID Arcview
	write(6,'("Result velocity output",/)')

!   
!	Escribimos resultados
    fname = 'Result_Vel.txt'
!
!	Normalizamos los resultados
	rmaxval = DMAX1(MAXVAL(acumuladoVel),1.d0)
	acumuladoVel = acumuladoVel / rmaxval
!
!	Escribimos la malla
	call GridOut(fname, acumuladoVel, mx, my, dx, xcorner, ycorner)
!
!
!
!	Calculamos la linea de máxima pendiente aprovechando la misma malla "Acumulado"
!
!	Inicializamos
	acumulado(:,:) = 0.d0
!
!	Diferenciamos el tipo de calculo
	if (itipo .EQ. 1) then
!
!		Llamamos a la fución que lo calcula y escribe el resultado en la malla "acumulado"
		Call MaxPath(ixini,iyini,acumulado,control,topo)
!
	else
!
!		Multiples sources
		do i = 1,mx
			do j = 1,my
!
!				Si se trata de u source calaculamos
				if(isources(i,j) .EQ. 1) then
					ixini = i
					iyini = j
!
!					Llamamos a la fución que lo calcula y escribe el resultado en la malla "acumulado"
					Call MaxPath(ixini,iyini,acumulado,control,topo)
!
				endif
!
			enddo
		enddo

	endif

!
!	Salida de resultados mediante GRID Arcview
	write(6,'("Stepest path output",/)')

!   
!	Escribimos resultados
	fname = 'MaxPend.txt'
!
!	Escribimos la malla
	call GridOut(fname, acumulado, mx, my, dx, xcorner, ycorner)
!
!
!
!	Salida de control para detección de sumideros mediante GRID Arcview
!
!
!	Salida de resultados mediante GRID Arcview
	write(6,'("Output of control matrix",/)')
!
!	No podemos escribir una malla de LOGICAL, pasamos a DOUBLE antes de excribir
	acumulado = merge(1.d0, 0.d0,control)
!   
!	Escribimos resultados
    fname = 'Control.txt'
!
!	Escribimos la malla
	call GridOut(fname, acumulado, mx, my, dx, xcorner, ycorner)
!
!
!	SIN ERROR
	goto 2000
!
!
!	ERROR
!
!	Se ha producido un error en la malla, no es igual a la de topo
1000 write(6,100)
100  format('Sourges grid and topo grid non corresponding',/)
!
!
!
2000 continue
!
!	Liberamos memoria
	DEALLOCATE(topo)				!Malla de topografia
	DEALLOCATE(topostat)			!Malla de pesos
	DEALLOCATE(acumulado)			!Malla de resultados
	DEALLOCATE(acumuladoVel)		!Malla de resultados de velocidad
	DEALLOCATE(control)				!Malla de control de borde y sumideros
	DEALLOCATE(isources)			!Malla de puntos de inicio
	DEALLOCATE(velocidades)			!Malla de velocidades
!
!
end program DebrisDice
!
!
!
