!****************************************************************************
!  Global data, parameters, and structures 
!****************************************************************************

module DebrisDiceGlobals
!
	implicit double precision (a-h,o-z)
!
!
!	Parameters of auxiliar data
	parameter (indexDzDx = 1)	!		#   aux(i,j,1) = dz/dx  (Value of Ix)
!
!
!   Global data
!
!
!   Enteros
	INTEGER :: nodata			!Nodata ESRI
	INTEGER :: itipo			!Coordenadas o malla
	INTEGER :: imaxpend		!Calculo de la maxima pendiente
	INTEGER :: iter			!Numero de iteraciones
	INTEGER :: iFricType		!Frictional(1) or viscous(2)
	INTEGER :: mx			!Numero nodos X
	INTEGER :: my			!Numero nodos Y
!
!
!   Booleanas
	LOGICAL :: Stopping
!
!   Reales
	REAL*8  :: xcorner		!Coordenada X esquina superior izq.
	REAL*8  :: ycorner		!Coordenada X esquina superior izq.
	REAL*8  :: dx			!Delta X
	REAL*8  :: dy			!Delta Y
	REAL*8  :: dd			!Delta diagonal
	REAL*8  :: rxini			!Xini
	REAL*8  :: ryini			!Yini
	REAL*8  :: abast			!Maxim abast
	REAL*8  :: amin			!Angulo minimo de parada
	REAL*8  :: rmu			!Coulomb friction factor
	REAL*8  :: rk			!Mass to drag ratio
	REAL*8  :: rm
	REAL*8  :: expm
	REAL*8  :: ctem
	REAL*8  :: cocm
	REAL*8  :: expBG		!Exponente
	REAL*8  :: viscBG		!Viscosidad, se da como dinamica (Pa·s)
	REAL*8  :: shBG			!ShearYield, damos la tension en Pa
	REAL*8  :: cte_BG
	REAL*8	:: densS		!Density
	REAL*8	:: bingDepth	!Aproximated depth used in viscous velocity calculation
!
!
!   Arrays
!
!   Enteros
	INTEGER,  DIMENSION(:,:),   ALLOCATABLE :: CoordCC          !Coordenadas de los nodos que son tipo 4 (CC)
	INTEGER  ,DIMENSION(5)                  :: method
!
!   Reales
	REAL*8,   DIMENSION(:,:,:), ALLOCATABLE :: q                !Matriz donde se almacenan las variables
	REAL*8,   DIMENSION(5)                  :: dtv
!
!
!   Texto
	CHARACTER*14 dummy
	CHARACTER*30 fname
!
!
end module DebrisDiceGlobals
