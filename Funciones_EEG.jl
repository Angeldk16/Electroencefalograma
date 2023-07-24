module Funciones_EEG

using Plots
import JLD.load , JLD.save
import Combinatorics.permutations
import Statistics.std , Statistics.cor , Statistics.mean
import CairoMakie.Point2f
import TopoPlots.eeg_topoplot
import CairoMakie.DataAspect
# Constantes
export OrdenMaximo
# Funciones para ambos
export VarsDeDoc , DatosDeTodos , Coordenadas_Y_Lados , COMS , ArrayParaHeatmap , PromVacios , CSDAZapfe , EliminaProhibidas , EliminaProhibidas , CentrosDeMasa , Dinamica_LV , IniciosYFinales , RellenaHuecos , ArreglaCoordenadas , DatosParaTopograma , DivideAnimacion
# Funciones para gráficas
export Registro_Completo , Animaciones_Trays
#En la parte derecha van las filas y en la izquierda las columnas
OrdenMaximo = Dict(
    "FP1" => [ 1 , 4 ],  "FPz" => [ 1 , 6 ], "FP2" => [ 1 , 8 ], 
    
    "AF7" => [ 2 , 4 ],  "AF3" => [ 2 , 5 ], "AFz" => [ 2 , 6 ], "AF4" => [ 2 , 7 ],  "AF8" => [ 2 , 8 ], 

    "F9" => [ 3 , 1 ],   "F7" => [ 3 , 2 ],  "F5" => [ 3 , 3 ],  "F3" => [ 3 , 4 ],   "F1" => [ 3 , 5 ], 
    "Fz" => [ 3 , 6 ],   "F2" => [ 3 , 7 ],  "F4" => [ 3 , 8 ],  "F6" => [ 3 , 9 ],   "F8" => [ 3 , 10 ], 
    "F10" => [ 3 , 11 ],  
    
    "FT9" => [ 4 , 1 ],  "FT7" => [ 4 , 2 ], "FC5" => [ 4 , 3 ], "FC3" => [ 4 , 4 ],  "FC1" => [ 4 , 5 ], 
    "FCz" => [ 4 , 6 ],
    "FC2" => [ 4 , 7 ],  "FC4" => [ 4 , 8 ], "FC6" => [ 4 , 9 ], "FT8" => [ 4 , 10 ], "FT10" => [ 4 , 11 ],
    
    "T7" => [ 5 , 2 ],   "C5" => [ 5 , 3 ],  "C3" => [ 5 , 4 ],  "C1" => [ 5 , 5 ],   "Cz" => [ 5 , 6 ], 
    "C2" => [ 5 , 7 ],   "C4" => [ 5 , 8 ],  "C6" => [ 5 , 9 ],  "T8" => [ 5 , 10 ], 
    
    "TP9" => [ 6 , 1 ],  "TP7" => [ 6 , 2 ], "CP5" => [ 6 , 3 ], "CP3" => [ 6 , 4 ],  "CP1" => [ 6 , 5 ], 
    "CPz" => [ 6 , 6 ], 
    "CP2" => [ 6 , 7 ],  "CP4" => [ 6 , 8 ], "CP6" => [ 6 , 9 ], "TP8" => [ 6 , 10 ], "TP10" => [ 6 , 11 ],
    
    "P9" => [ 7 , 1 ],   "P7" => [ 7 , 2 ],  "P5" => [ 7 , 3 ],  "P3" => [ 7 , 4 ],   "P1" => [ 7 , 5 ], 
    "Pz" => [ 7 , 6 ], 
    "P2" => [ 7 , 7 ],   "P4" => [ 7 , 8 ],  "P6" => [ 7 , 9 ],  "P8" => [ 7 , 10 ],  "P10" => [ 7 , 11 ],
    
    "PO7" => [ 8 , 4 ],  "PO3" => [ 8 , 5 ], "POz" => [ 8 , 6 ], "PO2" => [ 8 , 7 ],  "PO8" => [ 8 , 9 ],
    
    "O1" => [ 9 , 4 ],   "Oz" => [ 9 , 6 ],  "O2" => [ 9 , 8 ]
);
GaussianKernel=[0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067]
#El operador de Laplace-Lindenberg
LaplacianTerm1=[[0 1 0]; [1 -4 1]; [0 1 0]]
LaplacianTerm2=[[0.5 0 0.5]; [0 -2 0]; [0.5 0 0.5]]
LaplacianKernel=(1-1/3)*LaplacianTerm1+(1/3)*LaplacianTerm2

#************************ PRIMER SECCIÓN: FUNCIONES QUE USAN TANTO EL ELECTRO DE 64CH COMO EL DE 32CH  ***********************#

#------------------------------------------------ Obtiene variables ----------------------------------------------------------#
# Descripción: Obtiene las variables que se encuentran en los EEG de 64ch 
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                    Vars                                              Datos
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos               
function VarsDeDoc( Datos )
    Vars = Dict{ String , Any }()
    for i in 1:7
        Head = Datos[ i , : ][ 1 ]
        Head = split(Head,"[")[ 2 ]
        Head = String(split(Head,"]")[ 1 ])
        Dat = Datos[ i , : ][ 2 ]
        Vars[ Head ] = Dat
    end
    return Vars
end
#--------------------------------------------- Obtención de los datos del EEG ------------------------------------------------#
# Descripción: Esta función lo que hace es guardar todos los datos de cada uno de los canales en un 
#              diccionario con datos tipo Float32
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              Headers , Canales                          Datos , FilaHeaders , FilaDatos
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#                                        
function DatosDeTodos( Datos , FilaHeaders , FilaDatos )
    DatosH = Datos[ FilaHeaders , : ]
    Headers = []
    for i in 1:length( DatosH )
        temp = split( DatosH[ i ] , "]" )
        if length( temp[ 1 ]) > 1
            push!( Headers , temp[ 1 ] )
        end
    end
    Headers = String.(Headers)
    Canales = Dict()
    for j in 1:length( Headers )
        Canales[ Headers[ j ] ] = Float32.( Datos[ FilaDatos : end , j ] )
    end
    return Headers , Canales
end
#------------------------------------------------- Asigna Coordenadas --------------------------------------------------------#
# Descripción: Le asigna coordenadas a todos los canales y rellena los vacios con coordenadas prohibidas
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#       Coordenadas, lados, Prohibidas                          OrdenMaximo, Headers
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#                                        
function Coordenadas_Y_Lados( OrdenMaximo , Headers )
    All = collect( keys( OrdenMaximo ) )
    Posiciones = []
    for Þ in Headers
        if Þ in All
            push!( Posiciones , OrdenMaximo[ Þ ] )
        end
    end
    CutMax = Int.( zeros( 9 , 11 ) );
    for Pos in Posiciones
        CutMax[ Pos[ 1 ] , Pos[ 2 ] ] = 1
    end
    Dimensiones = size( CutMax )
    Acortar = [ ]
    for σ in 1 : Dimensiones[ 1 ]
        if length( unique( CutMax[ σ , : ] ) ) == 1 && CutMax[ σ , 1 ] != 1
            push!( Acortar , σ )
        end
    end
    if Acortar != Any[]
        CutMax2 = Int.( zeros( Dimensiones[ 2 ] ) )'
        if !(1 in Acortar )
            inicio = 0
            for ◉ in Acortar
                CutMax2 = vcat( CutMax2 , CutMax[ inicio + 1 : ◉ - 1 , : ] )
                inicio = ◉
            end
            if !( size( CutMax )[ 1 ] in Acortar )
                CutMax2 = vcat( CutMax2 , CutMax[ Acortar[ end ] + 1 : end , : ] )
            end
        else
            inicio = nothing
            for i in 2:length( Acortar )
                Δ = Acortar[ i ] - Acortar[ i - 1 ]
                if Δ > 1
                    inicio = Acortar[ i - 1 ]
                    break
                end
            end
            for ◉ in Acortar
                if ◉ < inicio
                    continue
                else
                    CutMax2 = vcat( CutMax2 , CutMax[ inicio + 1 : ◉ - 1 , : ] )
                    inicio = ◉
                end
            end
            if !( size( CutMax )[ 1 ] in Acortar )
                CutMax2 = vcat( CutMax2 , CutMax[ Acortar[ end ] + 1 : end , : ] )
            end
        end
        CutMax2 = CutMax2[ 2 : end , : ];
    else
        CutMax2 = CutMax
    end
    Acortar = [ ]
    for σ in 1 : Dimensiones[ 2 ]
        if length( unique( CutMax2[ : , σ ] ) ) == 1 && CutMax2[ 1 , σ ] != 1
            push!( Acortar , σ )
        end
    end
    if Acortar != Any[]
        CutMax3 = Int.( zeros( size( CutMax2 )[ 1 ] ) )
        if !(1 in Acortar )
            inicio = 0
            for ◉ in Acortar
                CutMax3 = hcat( CutMax3 , CutMax2[ : , inicio + 1 : ◉ - 1 ] )
                inicio = ◉
            end
            if !( size( CutMax2 )[ 2 ] in Acortar )
                CutMax3 = hcat( CutMax3 , CutMax2[ : , Acortar[ end ] + 1 : end ] )
            end
        else
            inicio = nothing
            for i in 2:length( Acortar )
                Δ = Acortar[ i ] - Acortar[ i - 1 ]
                if Δ > 1
                    inicio = Acortar[ i - 1 ]
                    break
                end
            end
            for ◉ in Acortar
                if ◉ < inicio
                    continue
                else
                    CutMax3 = hcat( CutMax3 , CutMax2[ : , inicio + 1 : ◉ - 1 ] )
                    inicio = ◉
                end
            end
            if !( size( CutMax2 )[ 2 ] in Acortar )
                CutMax3 = hcat( CutMax3 , CutMax2[ : , Acortar[ end ] + 1 : end ] )
            end
        end
        CutMax3 = CutMax3[ : , 2 : end ];
    else
        CutMax3 = CutMax2
    end
    CutMax = CutMax3
    lados = Int.( zeros( size( CutMax )[ 1 ] ) );
    for ⅁ in 1 : length( lados )
        ℑ = findall( x -> x == 1 , CutMax[ ⅁ , : ] )
        lados[ ⅁ ] = ℑ[ end ] - ℑ[ 1 ] + 1
    end
    NewValues = findall( x -> x == 1 , CutMax )
    NewCords = []
    for ∎ in 1 : length( NewValues )
        NewCord = [ NewValues[ ∎ ][ 1 ] , NewValues[ ∎ ][ 2 ] ]
        push!( NewCords , NewCord )
    end
    NewCords = sort( NewCords );
    Coordenadas = Dict()
    for ℜ in 1 : length( Posiciones )
        Canal = findfirst( Ϡ -> Ϡ == Posiciones[ ℜ ] , OrdenMaximo )
        if Canal in All
            Coordenadas[ Canal ] = NewCords[ ℜ ]
        end
    end
    ( Coordenadas , Prohibidas ) = Coordenadas_De_Vacios( CutMax , Coordenadas ) 
    return Coordenadas , lados, Prohibidas
end
#------------------------------------------------ Combinaciones de vecinos ---------------------------------------------------#
# Descripción: Nos da las combinaciones de +1 , 0 , 1 tales que sumadas a las coordenadas de un punto
#              dan como resultado sus vecinos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                 combinaciones                                  elementos
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#                Combinatorics   
function COMS( elementos )
    combinaciones = collect( permutations( elementos , 2 ) );
    extsup = [ 1 , 1 ]
    extinf = [ -1 , -1 ]
    push!( combinaciones , extsup )
    push!( combinaciones , extinf )
    return combinaciones
end
#------------------------------------------------ Creamos matriz de datos ----------------------------------------------------#
# Descripción: Toma los datos de los canales importantes y los arregla en forma de array para fines
#              de graficación y aplicación de procesos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              MVCuadrado                             lados , BIN01 , Prohibidas , Coordenadas
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#    
function ArrayParaHeatmap( lados , BIN01 , Prohibidas , Coordenadas )
    h = length( lados )
    l = maximum( lados )
    z = length( BIN01["O1"] )
    MVCuadrado = zeros( h , l , z )
    P = keys( Prohibidas )
    for k in 1:z
        for i in 1:h , j in 1:l
            ubi = [ i , j ]
            Can = findfirst( x -> x == ubi , Coordenadas )
            if !(Can in P)
                MVCuadrado[ i , j , k ] = BIN01[ Can ][ k ]
            end
        end
    end
    return MVCuadrado
end
#---------------------------------------------- Promediamos elementos vacios -------------------------------------------------#
# Descripción: Promedia los elementos que no tienen un valor
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              MVConPromedios                             Prohibidas , Coor , combs , Headers ,
#                                                                     BIN01 , MV
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#                Statistics                                        EncuentraVecinos
function PromVacios( Prohibidas , Coor , combs , Headers , BIN01 , MV )
    CoorPro = collect( keys( Prohibidas ) );
    VecinosDePro = Dict()
    for ν in 1:length( CoorPro )
        vecinos , indices = EncuentraVecinos( CoorPro[ ν ] , Coor , Prohibidas , combs , Headers )
        VecinosDePro[ CoorPro[ ν ] ] = vecinos
    end
    MVConPromedios = copy( MV )
    P = keys( Prohibidas )
    for k in 1:size( MV )[3]
        for i in 1:8 , j in 1:9
            ubi = [ i , j ]
            Can = findfirst( x -> x == ubi , Coor )
            if Can in P
                Vecinos = VecinosDePro[ Can ]
                Temp = []
                promed = 0
                for ξ in 1:length( Vecinos )
                    push!( Temp , BIN01[ Vecinos[ ξ ] ][ k ] )
                    promed = mean( Temp )
                end
                MVConPromedios[ i , j , k ] = promed
            end
        end
    end
    return MVConPromedios
end
#--------------------------------------------------- Obtenemos CSDA ----------------------------------------------------------#
# Descripción: Promedia los elementos que no tienen un valor
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              CSDA ( ∇ )                               Data , GaussianKernel , LaplacianKernel
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#                                                       GaussSuavizarTemporal , GaussianSmooth , 
#                                                            DiscreteLaplacian , UnNormGauss
"""
    CSDAZapfe( Data::Matrix{Float64} )
        -> GST::Matrix{Float16}, GS::Matrix{Float116}, CSD::Matrix{Float64}
        using GaussSuavizarTemporal, GaussianSmooth, DiscreteLaplacian
""" 
function CSDAZapfe( Data )
    #nChsl, nChsh, nFrs = size( Data );
    #side = Int( sqrt( nChs ) );
    #Data3D = reshape( Data, nChsl, nChsh, nFrs );
    ( μ, ν, ι )  = size( Data );
    # We apply a Temporal Gaussian smoothing ( this greatly affects the animations )
    Data3Plain = zeros( μ, ν, ι );
    for j = 1 : μ, l = 1 : ν
        channel = vec( Data[ j, l, : ] );
        Data3Plain[ j, l, : ] = GaussSuavizarTemporal( channel );
    end
    Φ = zeros( μ, ν, ι );
    ∇ = zeros( μ, ν, ι );
    # We spatially smooth the LFP with a two-dimensional Gaussian filter.
    # Later we obtain the dCSD.
    for τ = 1 : ι
        Φ[ :, :, τ ] = GaussianSmooth( Data3Plain[ :, :, τ ] );
        ∇[ :, :, τ ] = DiscreteLaplacian( Φ[ :, :, τ ] );
    end
    ∇ = -1 * ∇;
    return ∇
end  
#---------------------------------------------- Retiramos valores prohibidos -------------------------------------------------#
# Descripción: Elimina los valores de las coordenadas que no son canales
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              MVcsda                                     MVcsda2 , Prohibidas , Coordenadas
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#                                                    
function EliminaProhibidas( MVcsda2 , Prohibidas , Coordenadas )
    MVcsda = copy( MVcsda2 )
    P = keys( Prohibidas )
    for k in 1:size( MVcsda )[3]
        for i in 1:8 , j in 1:9
            ubi = [ i , j ]
            Can = findfirst( x -> x == ubi , Coordenadas )
            if Can in P
                MVcsda[ i , j , k ] = 0
            end
        end
    end
    return MVcsda
end
#----------------------------------------------- Obtenemos Centros de Masa ---------------------------------------------------#
# Descripción: Encuentra los elementos que superan algún umbral, encuentra las agrupaciones de 
#              dichos valores y calcula el centro de masa de cada agrupación
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              CMP , CMN                                  csda , factor
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#              Statistics                           ObtenComponentesyCM , ComponentesSP , vecindad8 ,
#                                                                     TiraOrillas
function CentrosDeMasa( csda , factor = 1 )
    ( ancho , alto , nmax ) = size( csda )
    scsd = std( csda )
    umbrsep = factor * scsd
    ( CMP, CMN ) = ObtenComponentesyCM( csda , 1 , nmax , umbrsep );
    return CMP , CMN
end
#------------------------------------------ Mapeamos el movimiento de los cm -------------------------------------------------#
# Descripción: Encuentra como se van moviendo los centros de masa anteriormente encontrados
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              TodasTrayectorias                             CMP , DistTol , TiempoTol
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#                                                      DaConsiderables , ECYB_Limites_Variables , 
#                                                     UltimasVelocidades  , dist2DVector , dist2D
function Dinamica_LV( Positivos2 , DistTol1 , TTol , Tmin )
    Positivos = copy(Positivos2)
    TotalFrames = length(Positivos)
    for k in 1:5
        Positivos[ TotalFrames + k ] = Array{Float64}(undef, 0, 3)
    end
    UltimoFrame = false
    pesomin=0
    TodasTrayectorias = Dict{Integer, Array{Any}}()
    ContadorTrayectorias = 1
    for tiempo in 1:TotalFrames
        if tiempo == TotalFrames
            UltimoFrame = true
        end
        ConjEnFrame = DaConsiderables( Positivos[ tiempo ], pesomin )
        NumConj = size( ConjEnFrame )[1]
        if NumConj > 0        
            for j in 1:NumConj
                if tiempo < TotalFrames
                    UltimoFrame = false
                end
                tiempoprimo = tiempo    
                HayMas = true
                ConjConFrame = [ transpose( ConjEnFrame[ j , : ] ) tiempo ]
                Trayectoria = ConjConFrame
                TExtra = 0
                while HayMas == true && TExtra <= TTol
                    if tiempoprimo >= TotalFrames
                        UltimoFrame = true
                    end                
                    if tiempoprimo < TotalFrames
                        ConjEnSigFrame = DaConsiderables( Positivos[ tiempoprimo + 1 ] 
                                                        , pesomin )
                        NumConjSig = size( ConjEnSigFrame )[ 1 ]
                    else
                        NumConjSig = 0
                        TExtra = TTol
                    end
                    if NumConjSig > 0
                        ( Trayectoria , Positivos , HayMas , TExtra , TodasTrayectorias 
                            , ContadorTrayectorias ) = ECYB_Limites_Variables( 
                            Trayectoria , Positivos , tiempoprimo , DistTol1 , HayMas, 
                            TExtra , TTol  , TodasTrayectorias , ContadorTrayectorias , Tmin )
                        if UltimoFrame == true
                            if size( Trayectoria )[1] > Tmin
                                TodasTrayectorias[ ContadorTrayectorias ] = Trayectoria
                                ContadorTrayectorias+=1
                                HayMas = false
                            end
                        end
                        tiempoprimo+=1
                    else
                        if TExtra == TTol
                            if size( Trayectoria )[ 1 ] > Tmin
                                TodasTrayectorias[ ContadorTrayectorias ] = Trayectoria
                                ContadorTrayectorias+=1
                                HayMas = false
                            else
                                HayMas = false
                            end
                        end
                        TExtra+=1
                        tiempoprimo+=1
                    end
                end # Cierra el while que ve si hay más 
            end # Cierra el for de cada conjunto
        end  # Cierra el if de si hay conjuntos en del frame
    end # Cierra el for de todos los frames
    return TodasTrayectorias
end
#------------------------------------------- Encuentra inicios y finales -----------------------------------------------------#
# Descripción: Crea dos diccionarios en donde se almacenan los tiempos en los que inicia y en los que 
#              termina cada trayectoria
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              Inicios , Finales                                       Trays
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#                                                     
function IniciosYFinales( Trays )
    Inicios = Dict{ Int , Float64 }()
    for i in 1:length( Trays )
        Inicios[ i ] = Trays[ i ][ 1 , end ]
    end
    Finales = Dict{ Int , Float64 }()
    for i in 1:length( Trays )
        Finales[ i ] = Trays[ i ][ end , end ]
    end
    return Inicios , Finales
end
#-------------------------------------------------- Rellena huecos -----------------------------------------------------------#
# Descripción: Cuando una trayectoria tiene muy pocos frames en los que no es continua los rellena 
#              con distancias y pesos intermedios
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                  TP                                         TP2 , Inicios , Finales
#      Las paqueterias que necesitamos         ///          Las funciones que necesitamos
#   
function RellenaHuecos( TP , Inicios , Finales )
    for i in 1:length( TP )
        if ( Finales[ i ] - Inicios[ i ] ) + 1 != size( TP[ i ] )[ 1 ]
            Frames = TP[ i ][ : , end ]
            for j in 1:( size( TP[ i ] )[ 1 ] - 1 )
                dif = Frames[ j + 1 ] - Frames[ j ]
                if dif > 1
                    Cabeza = TP[ i ][ 1 : j , : ]
                    Tronco = TP[ i ][ j + 1 : end , : ]
                    for k in 1:dif-1
                        h1 = (TP[ i ][ j , 1 ] + ( TP[ i ][ j+1 , 1 ] - TP[ i ][ j , 1 ] )*( k/dif ))
                        h2 = (TP[ i ][ j , 2 ] + ( TP[ i ][ j+1 , 2 ] - TP[ i ][ j , 2 ] )*( k/dif ))
                        h3 = (TP[ i ][ j , 3 ] + ( TP[ i ][ j+1 , 3 ] - TP[ i ][ j , 3 ] )*( k/dif ))
                        h4 = ( Inicios[ i ] + j + k - 1 )
                        Añadido = [ h1 h2 h3 h4 ]
                        Cabeza = vcat( Cabeza , Añadido )
                    end
                    TP[ i ] = vcat( Cabeza , Tronco )
                end
            end
        end
    end
    return TP
end
#--------------------------------------------- Transforma cuadrado a circulo -------------------------------------------------#
# Descripción: Las coordenadas cuadradas originales las transforma a circulares para que sean acordes 
#              al equipo con que se miden
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                    Trays                                         Trays , lados
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                                                          CuadACirc , polar2cartesian
function ArreglaCoordenadas( Trays , lados )
    for i in 1:length( Trays )
        # Primero recorremos los puntos en el eje x los cuales están en el rango [1,5]
        Trays[ i ][ : , 1 ] = Trays[ i ][ : , 1 ] .- ( ( length( lados ) - 1 ) / 2 + 1 )
        Trays[ i ][ : , 2 ] = Trays[ i ][ : , 2 ] .- ( ( maximum( lados ) - 1 ) / 2 + 1 )
        Trays[ i ][ : , 1 ] = Trays[ i ][ : , 1 ] ./ length( lados ) .* maximum( lados )
        Trays[ i ][ : , 1 : 2 ] = CuadACirc( Trays[ i ][ : , 1 : 2 ] )
    end
    return Trays
end
#-------------------------------------------- Cambia el formato de los datos -------------------------------------------------#
# Descripción: Cambia el formato de las posiciones de los canales a tipo Point2f  y crea diccionarios necesarios para graficar
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#            pesos , CT , CoorCSDA                      Coordenadas , Prohibidas , MVcsdaListo
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                  CairoMakie                                    ArrgelaCoordenadas
function DatosParaTopograma( Coordenadas , Prohibidas , MVcsdaListo , lados )
    CT = Dict( 1 => Float64.(getindex.( collect( values( Coordenadas ) ), [ 1 2 ] ) ) );
    PT = collect( values( Prohibidas ) );
    CTT = zeros( size( CT[1] )[ 1 ] - length( PT ) , 2 );
    τ = 0 
    for i in 1:size(CT[1])[1]
        ubi = CT[1][i,:]
        if !(ubi in PT)
            τ = τ + 1
            CTT[ τ , : ] = ubi
        end
    end
    CT = Dict( 1 => CTT )
    pesos = []
    for k in 1:size(MVcsdaListo)[3]
        Is = MVcsdaListo[:,:,k]
        W = []
        for i in 1:size(CT[1])[1]
            ubi = [ Int.(CT[1][i,1]) , Int.(CT[1][i,2]) ]
            if !(ubi in PT)
                w = Is[ Int.(CT[1][i,1]) , Int.(CT[1][i,2]) ]
                push!( W , w )
            end
        end
        W = Float64.( W )
        push!( pesos , W )
    end
    CT = ArreglaCoordenadas( CT , lados )
    CoorCSDA = []
    for i in 1:size(CT[1])[1]
        Ξ = [ CT[1][ i , 1 ] , CT[1][ i , 2 ] ]
        push!( CoorCSDA , Ξ )
    end
    CoorCSDA = Point2f.( CoorCSDA )
    return pesos , CT , CoorCSDA
end
#------------------------------------------------ Selecciona bines de tiempo -------------------------------------------------#
# Descripción: Debido a que algunas animaciones serán algo pesados no es recomendable hacer una sola animación muy pesada por   # lo que esta función divide los frames en partes iguales y toma ciertos bines al azar
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#               cachos , Bines                                      PC , pesos
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                                                   
function DivideAnimacion( PC , pesos )
    cachos = 0
    finish = false
    PCI = deepcopy( PC )
    while finish == false
        TT = length( pesos )
        if TT%PC == 0
            cachos = PC
            finish = true
        end
        if PC > PCI * 2
            if floor(TT/PC)*PC < length( pesos )
                cachos = PC
                finish = true
            end
        end
        PC+=1
    end
    Bines = rand( 1 : cachos , 5 )
    Bines = unique( Bines )
    return cachos , Bines
end

#**************************** SEGUNDA SECCIÓN: FUNCIONES AUXULIARES PARA LAS FUNCIONES PRINCIPALES ***************************#

#------------------------------------------- Complemento de: Coordenadas_Y_Lados ---------------------------------------------#
# Descripción: Los segmentos del array que no pertenecen a un canal les asigna coordenadas y los nombra
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#          Coordenadas, Prohibidas                             CutMax , Coordenadas
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#                                        
function Coordenadas_De_Vacios( CutMax , Coordenadas )  
    Prohibidas = Dict()
    ( l , a ) = ( size( CutMax )[ 1 ] , size( CutMax )[ 2 ] )
    PRN = "VACIO"
    ND = length( string( length( findall( x -> x == 0 , CutMax ) ) ) ) + 1
    ᶲ = 1
    for Ŋ in 1 : l
        for Ƶ in 1 : a
            Casilla = CutMax[ Ŋ , Ƶ ]
            if Casilla == 0
                Name = PRN * string( lpad( ᶲ , ND , "0" ) )
                Coordenadas[ Name ] = [ Ŋ , Ƶ ]
                Prohibidas[ Name ] = [ Ŋ , Ƶ ]
                ᶲ = ᶲ + 1
            end
        end
    end
    return Coordenadas , Prohibidas         
end
#----------------------------------------------- Complemento de: PromVacios --------------------------------------------------#
# Descripción: Encuentra la 8-vecindad de un elemento y dice sus coordenadas 
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#           Vecinos , indices                      canal , Coordenadas , Prohibidas , combinaciones , 
#                                                                      Headers
function EncuentraVecinos( canal , Coordenadas , Prohibidas , combinaciones , Headers )
    ValCor = values( Coordenadas )
    ValPro = values( Prohibidas )
    CoorVecinos = []
    CoorCanal = Coordenadas[ canal ]
    for tup in combinaciones
        temp = CoorCanal + tup
        if temp in ValCor && !(temp in ValPro)
            push!( CoorVecinos , temp )
        end
    end
    Vecinos = findall( x -> x in CoorVecinos , Coordenadas)
    indices = []
    for canal in Vecinos
        name = string(canal)
        index = findfirst( x -> x == name, Headers )
        push!( indices , index )
    end
    return Vecinos , indices
end
#----------------------------------------------- Complemento de: CSDAZapfe ---------------------------------------------------#
# Descripción: Un suavizado Gaussiano temporal. Esto es escencialmente un filtro pasabajos. Depende 
#              implicitamente de la frecuencia de muestreo. sigma esta medido en pixeles, es la 
#              desviacion estandar de nuestro kernel. El medioancho de nuestra ventana seran 3*sigma
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#               result                                             Datos , Sigma 
function GaussSuavizarTemporal( Datos , Sigma = 3 )
    medioancho = ceil( Sigma * 3 )
    colchon = ones( medioancho )
    result = zeros( size( Datos ) )
    datoscolchon = vcat( colchon * Datos[ 1 ] , Datos , colchon * Datos[ end ] )
    kernel = map( x -> UnNormGauss( x , Sigma ), collect( -medioancho : medioancho ) )
    kernel = kernel/( sum( kernel ) )
    for t = medioancho + 1 : length( Datos ) + medioancho
        result[ t - medioancho ] = sum( datoscolchon[ t - medioancho : t + medioancho ] .* kernel )
    end
    return result
end
#----------------------------------------------- Complemento de: CSDAZapfe ---------------------------------------------------#
# Descripción: Realiza un suavizado gaussiano espacial
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                  result                                                Datos
function GaussianSmooth( Datos )
    tamanodatos = size( Datos )
    result = zeros( tamanodatos )
    temp = copy( Datos )
    ( mu , lu ) = size( Datos )
    #Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    ## Okey, parece que los imbeciles de rioarriba cambiaron la sintaxis de
    # rebanadas de matriz. Ahora CUALQUIER rebanada de matriz es colvec.
    arriba = reshape( temp[ 1 , : ] , ( 1 , lu ) )
    abajo = reshape( temp[ end , : ] , ( 1 , lu ) )
    arr3 = vcat( arriba , arriba , arriba )
    aba3 = vcat( abajo , abajo , abajo )   
    temp = vcat( arr3 , temp, aba3 ) 
    for j=1:3
        temp = hcat( temp[ : , 1 ] , temp , temp[ : , end ] )
    end
    for j = 4 : tamanodatos[ 1 ] + 3 , k = 4 : tamanodatos[ 2 ] + 3
        aux = temp[ j - 3 : j + 3 , k - 3 : k + 3 ]
        result[ j - 3 , k - 3 ] = sum( GaussianKernel .* aux )
    end
    #Esta convolución no respeta norma L2
    #result=result*maximum(abs(Datos))/maximum(abs(result))
    return result
end
#---------------------------------------------- Complemento de: CSDAZapfe ----------------------------------------------------#
# Descripción: Aplica el operador Laplaciano a la matriz
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                  result                                                Datos
function DiscreteLaplacian( Datos )
    temp = copy( Datos )
    ( mu , lu ) = size( Datos )
    izq = reshape( temp[ 1 , : ] , ( 1 , lu ) )
    der = reshape( temp[ end , : ] , ( 1 , lu ) ) 
    #Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    temp = vcat( izq , temp , der )
    temp = hcat( temp[ : , 1 ] , temp , temp[ : , end ] )
    largo , ancho = size( temp )
    aux = Array{ Float32 }( undef , 3 , 3 )
    result=zeros( size( temp ) )    
    for j = 2 : largo - 1 , k = 2 : ancho - 1
        #los indices van primero, "renglones", luego "columnas", etc
        aux = temp[ j - 1 : j + 1 , k - 1 : k + 1 ]
        result[ j , k ] = sum( LaplacianKernel .* aux )
    end
    #DO  Crop the borders
    result = result[ 2 : end - 1 , 2 : end - 1 ]
    return result
end
#---------------------------------------------- Complemento de: CSDAZapfe ----------------------------------------------------#
function UnNormGauss( x , sigma )
    return exp( - x * x / ( 2 * sigma ) )
end
#-------------------------------------------- Complemento de: CentrosDeMasa --------------------------------------------------#
# Descripción: Obtiene los componentes de todos los conjuntos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#           CMPositivo, CMNegativo                        Datos , tini , tfini , epsilon
function ObtenComponentesyCM( Datos::Array , tini = 1 , tfini = tmax , epsilon = 1.0 )
    #CSD ahora no tiene orillas. Asi que toca adaptarse.
    ( alto , ancho , lu ) = size( Datos )
    #la cantidad minima de pixeles que tiene que tener un componente para
    #que lo tomemeos en cuenta
    tamano = 3
    CMPositivo = Dict{ Int, Array }()
    CMNegativo = Dict{ Int, Array }()
    for t = tini : tfini
        ActividadNegativa = Array{ Int16 }[]
        ActividadPositiva = Array{ Int16 }[]
        SpikeCountPositivo = zeros( alto , ancho )
        SpikeCountNegativo = zeros( alto , ancho )
        #Separamos pixeles positivos y negativos
        for j = 1 : alto , k = 1 : ancho
            if( Datos[ j , k , t ] < -epsilon )
                push!( ActividadNegativa , [ j , k ] )
                SpikeCountNegativo[ j , k ]+=1
            elseif( Datos[ j , k , t ] > epsilon )
                push!( ActividadPositiva , [ j , k ] )
                SpikeCountPositivo[j,k]+=1
            end
        end
            #Primero Negativo
        componentesneg = ComponentesSP( ActividadNegativa )
        centrosdemasaneg=[[0 0 0];]
        for p in componentesneg
            mu = length( p )        
            if mu > tamano
                masa = 0.00
                x = 0.00
                y = 0.00
                for q in p
                    j = q[ 1 ]
                    k = q[ 2 ]
                    masalocal = Datos[ j , k , t ]
                    masa+=masalocal
                    x+=k*masalocal          
                    y+=j*masalocal
                end
                x/=masa            
                y/=masa
                A=[x y masa]               
                centrosdemasaneg = vcat( centrosdemasaneg , A )
            end
        end
        centrosdemasaneg = centrosdemasaneg[ 2 : end , : ]
        CMNegativo[ t ] = centrosdemasaneg
        ##### Ahora lo posittivo (fuentes)
        componentespos = ComponentesSP( ActividadPositiva )               
        centrosdemasapos=[[0 0 0];]
        for p in componentespos
            mu = length( p )
            if mu > tamano
                masa = 0.00
                x = 0.00
                y = 0.00
                for q in p
                    j = q[ 1 ]
                    k = q[ 2 ]
                    masalocal = Datos[ j , k , t ]
                    masa+=masalocal
                    x+=k*masalocal
                    y+=j*masalocal
                end
                x/=masa 
                y/=masa
                A=[x y masa]
                centrosdemasapos = vcat( centrosdemasapos , A )
            end
        end
        centrosdemasapos = centrosdemasapos[ 2 : end , : ]       
        CMPositivo[ t ] = centrosdemasapos
    end
    return ( CMPositivo , CMNegativo )
end
#----------------------------------------------- Complemento de: CentrosDeMasa -----------------------------------------------#
# Descripción: Obtiene los componentes que superan un umbral
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                componentes                                        DatosSignados  
function ComponentesSP( DatosSignados::Array )
    #Single pass method for Disjoint Components.
    lista = copy( DatosSignados )
    componentes = Set{ Any }()
    while( length( lista ) != 0 )
        x = pop!(lista) #arranca el ULTIMO elemento de la lista
        listaprofundeza = Array{ Int64 }[]
        componentecurlab = Array{ Int64 }[]
        push!( listaprofundeza , x ) #Pone elementos al FINAL de la lista
        push!( componentecurlab , x )    
        profundidad = 0
        while ( ( length( listaprofundeza ) != 0 ) && profundidad < 1000)
            y = pop!( listaprofundeza )
            for v in vecindad8( y )
                if in( v , lista ) 
                    deleteat!( lista , indexin( Any[ v ] , lista ) )
                    push!( listaprofundeza , v ) 
                    profundidad+=1
                    push!( componentecurlab , v )
                end
            end
        end
        push!( componentes , componentecurlab )    
    end
    return componentes
end
#---------------------------------------------- Complemento de: CentrosDeMasa ------------------------------------------------#
# Descripción: Encuentra la 8-vecindad
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#           8-vecindad ( result )                                       punto
function vecindad8( punto::Array )
    j = punto[ 1 ]
    k = punto[ 2 ]
    vecindad = Set{ Array{ Int64 , 1 } }()
    push!( vecindad , [ j - 1 , k - 1 ] )
    push!( vecindad , [ j - 1 , k ] )
    push!( vecindad , [ j - 1 , k + 1 ] )
    push!( vecindad , [ j , k - 1 ] )
    push!( vecindad , [ j , k + 1 ] )
    push!( vecindad , [ j + 1 , k - 1 ] )
    push!( vecindad , [ j + 1 , k ] )
    push!( vecindad , [ j + 1 , k + 1 ] )
    vecindad = TiraOrillas( vecindad )
    return vecindad
end
#--------------------------------------------- Complemento de: CentrosDeMasa -------------------------------------------------#
# Descripción: Quita las orillas de una matriz cuadrada
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                  result                                             Puntos
function TiraOrillas( Puntos::Set )
    result=Set([])
    for p in Puntos
        if !(p[1]==0 || p[2]==0 || p[1]==65 ||  p[2]==65)
            push!( result , p )
        end
    end
    return result
end
#--------------------------------------------- Complemento de: Dinamica_LV ---------------------------------------------------#
# Descripción: Nos dará los conjuntos en el frame que superen el peso mínimo
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#             Considerables                                     Frame , PesoMin
function DaConsiderables( Frame , PesoMin )
    NF = size( Frame )[ 1 ]
    if NF > 0
        Pesos = Frame[ : , 3 ]
        EncuentraConsiderables = findall( abs.( Pesos ) .> PesoMin )
        Considerables = Frame[ EncuentraConsiderables , : ]
    else
        Considerables = Frame
    end
    return Considerables
end
#--------------------------------------------- Complemento de: Dinamica_LV ---------------------------------------------------#
# Descripción: Encuentra la continuación de la trayectoria, la añade y la borra de los datos
function ECYB_Limites_Variables( Trayectoria , Positivos , tiempo , DistTol1 , 
    HayMas , TExtra , TTol , TodasTrayectorias , ContadorTrayectorias , Tmin )
    if size(Trayectoria)[1] < 5
        DistTol = DistTol1
    else
        Locs = Trayectoria[ ( end - 4 ) : ( end ) , [ 1 , 2 ] ]
        Ts = Trayectoria[ ( end - 4 ) : ( end ) , 4 ]
        ArrVel = []
        for frame in 2:5
            ( V , ΔT ) = UltimasVelocidades( Ts , Locs , frame  )
            push!( ArrVel , V )
        end
        Ts = Ts[2:end]
        X2 = zeros( 4 ,2);
        X2[:,1] = transpose(Ts);  
        X2[:,2] .= 1.0
        ArrVel = Float64.( ArrVel )
        coeff_pred = X2\ArrVel
        DistTol = ( coeff_pred[1]*(Ts[end]+1)+coeff_pred[2])*1.5
        if DistTol < DistTol1/5
            DistTol = DistTol1/5
        end
    end
    ultimo = Trayectoria[end,:]
    distancias = dist2DVector(ultimo,Positivos[ tiempo + 1 ])
    MasCercano = argmin(distancias)
    DistMasCercano = minimum(distancias)
    if DistMasCercano < DistTol
        Temporal = [ transpose(Positivos[ tiempo + 1 ][ MasCercano , : ]) tiempo + 1 ]
        Cadena = vcat( Trayectoria , Temporal)
        Positivos[ tiempo + 1 ] = Positivos[ tiempo + 1 ][ 1:end .!= MasCercano, : ]
        TExtra = 0
    else
        Cadena = Trayectoria
        if TExtra == TTol
            HayMas = false
            if size( Trayectoria )[ 1 ] > Tmin
                TodasTrayectorias[ ContadorTrayectorias ] = Trayectoria
                ContadorTrayectorias+=1
            end
        end
        TExtra+=1
    end
    return Cadena , Positivos , HayMas , TExtra , TodasTrayectorias , ContadorTrayectorias
end
#--------------------------------------------- Complemento de: Dinamica_LV ---------------------------------------------------#
# Descripción: Encuentra las ultimas velocidades de la trayectoria 
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                   V , ΔT                                     Ts , Locs , frame
function UltimasVelocidades( Ts , Locs , frame  )
    x = Locs[ frame - 1 , : ]
    y = Locs[ frame , : ]
    ΔLoc = dist2D(x,y)
    ΔT = Ts[ frame ] - Ts[ frame  - 1 ]
    V = ΔLoc / ΔT
    return V , ΔT
end
#--------------------------------------------- Complemento de: Dinamica_LV ---------------------------------------------------#
# Descripción: encuentra las distancias entre las posiciones x y y de un vector con todos los
#              elementos de una matriz 
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                   dist                                            v , M
# v es el vector donde v[1] y v[2] son las posiciones en x y y, M es la matriz
function dist2DVector(v,M)
    dist = sqrt.( ( v[ 1 ] .- M[ :, 1 ] ) .^ 2 .+ ( v[ 2 ] .- M[ :, 2 ] ) .^ 2 );
    return dist
end
#-------------------------------------------- Complemento de: Dinamica_LV ----------------------------------------------------#
# Descripción: encuentra las distancias entre 2 puntos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                   dist                                            v , M
function dist2D( x , y )
    result = sqrt( ( x[ 1 ] - y[ 1 ] )^2 + ( x[ 2 ] - y[ 2 ] )^2 )
    return result
end
#--------------------------------------- Complemento de: ArreglaCoordenadas --------------------------------------------------#
# Descripción: Transforma cuadrados a circulos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#              PuntosPolares                                           puntos
function CuadACirc( puntos )
    PuntosPolares = zeros( size( puntos )[ 1 ] , 2 )
    for i in 1:size( puntos )[ 1 ]
        if abs(puntos[i,1]) > abs(puntos[i,2]) # Si el ancho es más grande que el largo
            r = abs(puntos[i,1])
            signoA = puntos[i,1]/abs(puntos[i,1])# positivo a la derecha, sino a la izquierda
            signoL = puntos[i,2]/abs(puntos[i,2])#positivo arriba, sino abajo
            θ = puntos[i,2]/r*(π/4)  #transforma las coordenadas del punto menor a radianes
            if signoA < 0
                θ = θ + (π)
            end
        elseif abs(puntos[i,2]) > abs(puntos[i,1])
            r = abs(puntos[i,2])
            signoA = puntos[i,1]/abs(puntos[i,1])#positivo a la derecha, sino a la izquierda
            signoL = puntos[i,2]/abs(puntos[i,2])#positivo arriba, sino abajo
            θ = puntos[i,1]/r*(π/4) #transforma las coordenadas del punto menor a radianes
            if signoL > 0 
                θ = -θ + (π/2)
            elseif signoL < 0 
                θ = -θ + (3*π/2)
            end
        else
            r = abs(puntos[i,2])
            if r > 0
                signoA = puntos[i,1]/abs(puntos[i,1])
                signoL = puntos[i,2]/abs(puntos[i,2]) #positivo arriba, sinp abajo
                if signoL > 0 
                    if signoA < 0
                        θ = puntos[i,2]/r*(π/4)+(π/2) # Extremo superior izquierdo
                    else
                        θ = puntos[i,2]/r*(π/4) # Extremo superior derecho
                    end
                else
                    if signoA < 0
                        θ = puntos[i,2]/r*(π/4)+(6*π/4) # Extremo inferior izquierdo
                    else
                        θ = puntos[i,1]/r*(π/4)+(6*π/4) # Extremo inferior derecho
                    end
                end
            else
                θ = 0
            end
        end
        ( x , y ) = polar2cartesian(r, θ)
        PuntosPolares[ i , 1 ] = x
        PuntosPolares[ i , 2 ] = y
    end
    return PuntosPolares
end
#------------------------------------------ Complemento de: ArreglaCoordenadas -----------------------------------------------#
# Descripción: Transforma coordenadas polares a cartesianas
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                 x, y                                                  r, θ
function polar2cartesian( r , θ )
    x = r * cos( θ )
    y = r * sin( θ )
    return ( x, y )
end

#******************************************* TERCER SECCIÓN: FUNCIONES PARA GRAFICAR *****************************************#

#--------------------------------------------------- Imagen del registro -----------------------------------------------------#
# Descripción: Grafica los datos del voltaje en todo el registro del experimento
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#            Figura del registro                         Headers, Canales, ruta_resultados
#      Las paqueterias que necesitamos      ///      Las funciones que necesitamos
#                   Plots                    
function Registro_Completo( Headers , Canales , ruta_resultados )
    A = nothing
    for i in 1 : length( Headers )
        if i == 1
            A = plot( Canales[ Headers[ i ] ] , label = Headers[ i ] , legend_position =:right , 
                wsize = (1000,500))
        else
            A = plot!( Canales[ Headers[ i ] ] .- 100 * ( i - 1 ) , label = Headers[ i ] , 
                legend_position =:right , wsize = (1000,500))
        end
    end
    Records_Image = joinpath( ruta_resultados , "Registros.png" );
    savefig( A , Records_Image );
    Listo = println( "La imagen está lista" )
    return Listo 
end
#------------------------------------------------ Crea Gif's de trays y csda -------------------------------------------------#
# Descripción: Crea un gif donde la parte izquierda muestra las trayectorias de los centros de masa y la parte derecha muestra un heatmap del csda
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                                                                      muchas 
#      Las paqueterias que necesitamos       ///           Las funciones que necesitamos
#            Plots , Statistics                 GeneraColores , TraysFrameXFrame , DinamicaEnFrame ,
#                                                                  GraficaDinamica
function Animaciones_Trays(extra ,cachos , nd , longtray , fr , MVcsda , ruta_resultados , IniciosP , FinalesP , TraysPosP , 
         CT , IniciosN , FinalesN , TraysNegP , CTPs , CTNs , bines, Headers = nothing, Canales = nothing , pesos = nothing ,
         CoorCSDA = nothing )
    CTP = GeneraColores( CTPs , TraysPosP );
    CTN = GeneraColores( CTNs , TraysNegP );
    minx , maxx = (minimum( CT[ 1 ][ : , 1 ] ) - extra, maximum( CT[ 1 ][ : , 1 ] ) + extra)
    miny , maxy = (minimum( CT[ 1 ][ : , 2 ] ) - extra, maximum( CT[ 1 ][ : , 1 ] ) + extra)
    tamaño = size( MVcsda )[ 3 ] / cachos
    nomcsda = "Analisis_"
    xylims = [ minx maxx miny maxy ]
    if maxx > abs( minx )
        minx = - maxx
    else
        maxx = - minx
    end
    if maxy > abs( miny )
        miny = - maxy
    else
        maxy = - miny
    end
    ruta_gifs = joinpath( ruta_resultados * "_GIFs" );
    mkpath( ruta_gifs );
    cd( ruta_gifs )
    # Sección necesaria para los EEG
    if Headers != nothing
        NameEEG = "EEG_TOPO"    
        tamaño = floor( length( pesos ) / cachos )
        N0s = Int( length( string( Int( tamaño ) ) ) )
        NGs = length( string( cachos ) )
        Cords = CoorCSDA
        Axs = (aspect=DataAspect(),)
        Cont = (color=:black, linewidth=0.1) 
    end
    for i in bines
        ruta_pngs = nothing
        if Headers != nothing
            ruta_pngs = joinpath( ruta_resultados , "Temporal" )
            mkpath( ruta_pngs )
            start = Int( ( tamaño * ( i - 1) + 1 ) )
            finish = Int( ( tamaño * i ) )
            for ν in 1:tamaño
                NombreFig = string( NameEEG , lpad( Int( ν ) , N0s , "0" ) , ".png" )
                DirToSave = joinpath( ruta_pngs , NombreFig )
                x = eeg_topoplot(pesos[Int(start + ν - 1)]; positions=Cords, contours=Cont, axis=Axs)
                save( DirToSave , x ; px_per_unit = 2 )
            end
        end
        start = Int(( tamaño * ( i - 1) + 1 ))
        finish = Int(( tamaño * i ))
        S = std( Float64.( MVcsda[:,:, start:finish ] ) )
        limites = ( -4*S , 4*S )
        ( ImXGifTrays , ImXGifCSDA , ImXGifRecord , ImXGifTopos ) = TraysFrameXFrame( longtray , start , finish , IniciosP , 
                                       FinalesP , TraysPosP , CTP , xylims, IniciosN , FinalesN ,
                                       TraysNegP , CTN , fr , MVcsda , limites, Headers , Canales ,
                                       ruta_pngs )
        anima = @animate for i in 1:length( ImXGifCSDA )
            time = string( ( start + i ) / fr)
            if Headers != nothing
                C = plot(ImXGifTopos[ i ] , ImXGifRecord[ i ] , ImXGifTrays[ i ] , ImXGifCSDA[ i ] , wsize =(800 , 800) , 
                         Layout = (2,2) , title = time )
            else
                C = plot(ImXGifTrays[ i ] , ImXGifCSDA[ i ] , wsize =(800 , 400) , Layout = (1,2) , title = time )
            end
        end
        NomGifCSDA = string( nomcsda , lpad( i , nd , "0") , ".gif" )
        NomGifTodo = string( "Analisis_Completo_" , lpad( i , nd , "0") , ".gif" )
        if Headers != nothing
            gif(anima, NomGifTodo , fps = 60 );
        else
            gif(anima, NomGifCSDA , fps = 60 );
        end
        println("Se guardó el registro ",i," de ",cachos)
    end
    Listo = "El proceso terminó"
    return Listo
end
#------------------------------------------------- Crea series de colores ----------------------------------------------------#
# Descripción: Para que en las gráficas no se confundan las diferentes gráficas tenemos que generar distintos colores para luego asignarlos
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                    CTP                                             CTs , TP2
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
# 
function GeneraColores( CTs , TP2 )
    CTP = Dict()
    LG = length(CTs)
    cc = 0
    for i in 1:length(TP2)
        cc+=1
        CTP[ i ] = CTs[cc]
        if cc == LG
            cc = 0
        end
    end
    return CTP
end
#--------------------------------------------- Añade la gráfica de una trayectoria -------------------------------------------#
# Descripción: Añade la trayectoria actual de un conjunto a la grafica del frame
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                    A                                                muchas
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                  Plots                        
function GraficaDinamica( c , frame , longtray , Origen , inicio , final , cg , xylims , fr)
    if frame - inicio + 1  < longtray
        x = Origen[ 1 : Int( frame - inicio + 1 ) , 1 ]
        y = Origen[ 1 : Int( frame - inicio + 1 ) , 2 ]
    elseif frame < final - longtray && frame >= longtray
        x = Origen[ Int( frame - inicio + 1 ) : Int( frame - inicio + 1 + longtray ) , 1 ]
        y = Origen[ Int( frame - inicio + 1 ) : Int( frame - inicio + 1 + longtray ) , 2 ]
    else
        x = Origen[ Int( frame - inicio + 1 ) : end , 1 ]
        y = Origen[ Int( frame - inicio + 1 ) : end , 2 ]
    end
    if c == 1
        A = plot!(
            x , y , color = cg
            ,line=(:dot,1), lab=frame/fr,
            marker=([:square :d],5,0.8,stroke(3,:gray)),  
            xlims = ( xylims[1] , xylims[2] ) , ylims = ( xylims[3] , xylims[4] ),
            leg=false
            )
    else 
        A = plot!(
            x , y , color = cg
            ,line=(:dot,1), lab=frame/fr,
            marker=([:square :d],5,0.8,stroke(3,:gray)),  
            xlims = ( xylims[1] , xylims[2] ) , ylims = ( xylims[3] , xylims[4] ),
            leg=false
            )
    end
    return A
end
#--------------------------------------------- Crea heatmap de pocos canales -------------------------------------------------#
# Descripción: Función de ayuda para mostrar canales que superan cierto valor
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                   Map                                 MVCuadradoCompleto , frame , umbral
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                  Plots                                  
function MapaSinCeros( MVCuadradoCompleto , frame , umbral )
    MVCuadradoCompleto2 = copy(MVCuadradoCompleto)
    Dimensiones = size( MVCuadradoCompleto2[ : , : , frame ] )
    for i in 1:Dimensiones[1], j in 1:Dimensiones[2]
        if MVCuadradoCompleto2[ i , j , frame ] < umbral
            MVCuadradoCompleto2[ i , j , frame ] = umbral
        end
    end
    Map = heatmap( MVCuadradoCompleto2[ : , : , frame ] )
    return Map
end
#------------------------------------------- Crea imagenes de todos los frames -----------------------------------------------#
# Descripción: Crea dos arrays con las imagenes del CSDA y de las trayectorias en cada frame
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#          ImXGifTrays , ImXGifCSDA                                   muchas
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                                                                 DinamicaEnFrame
function TraysFrameXFrame( longtray , F1 , FF , IniciosP , FinalesP , TP2 , CTP , xylims , IniciosN , FinalesN , TN2 , CTN ,
                           frecuencia , MVcsda , limites, Headers, Canales , ruta_pngs )
    ActivosP = Dict()
    ActivosN = Dict()
    PON = [ :square :triangle ]
    ImXGifTrays = []
    ImXGifCSDA = []
    ImXGifRecord = []
    ImXGifTopos = []
    for frame in F1:FF
        ( ActivosP, ActivosN, A, B, C, D ) = DinamicaEnFrame( frame , IniciosP , IniciosN , FinalesP , FinalesN , ActivosP ,                                   ActivosN , TP2 , TN2 , longtray , CTP , CTN , xylims , frecuencia , MVcsda , limites, 
                                Headers , Canales , ruta_pngs , F1, FF)
        push!( ImXGifTrays , A )
        push!( ImXGifCSDA , B )
        if Headers != nothing
            push!( ImXGifRecord , C )
            push!( ImXGifTopos , D )
        end
    end
    return ImXGifTrays , ImXGifCSDA , ImXGifRecord, ImXGifTopos
end
#------------------------------------------ Grafica las trayectorias en el frame ---------------------------------------------#
# Descripción: Crea la gráfica donde se muestran todas las trayectorias que están sucediendo en ese frame
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#        ActivosP , ActivosN , A , B                                    muchas
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                 Plots                                           GraficaDinamica
function DinamicaEnFrame( frame , IniciosP , IniciosN , FinalesP , FinalesN , ActivosP , ActivosN ,TP2 , TN2 , longtray , 
                          CTP , CTN , xylims , frecuencia , MVcsda , limites, Headers , Canales , ruta_pngs, F1, FF
                          )
    TempInP = findall(x->x==frame,IniciosP)
    for i in TempInP
        ActivosP[ i ] = [ TP2[i][:,1] , TP2[i][:,2] ]
    end
    TempFinP = findall(x->x==frame,FinalesP)
    for i in TempFinP
        delete!(ActivosP,i)
    end
    ValActP = sort(collect(keys(ActivosP)));
    ValActN = sort(collect(keys(ActivosN)));
    radio = maximum( abs.( xylims ) )
    PuntosX = -1:0.01:1
    xs = [ ]
    ys = [ ]
    for eq in PuntosX
        h = sqrt( 1 - eq^2 )
        push!( xs , eq )
        push!( ys , h )
    end
    xs = vcat( xs , reverse( xs ) )
    ys = vcat( ys , -ys )
    A = plot( radio*xs[ 1 : end ] , radio*ys[ 1 : end ] , leg=false , color =:black 
              , xlims = ( xylims[1] , xylims[2] ) , ylims = ( xylims[3] , xylims[4] ) )
    c = 0
    for i in 1:length(ValActP)
        c+=1
        A = GraficaDinamica(c,frame,longtray,TP2[ValActP[i]]
            ,IniciosP[ValActP[i]],FinalesP[ValActP[i]],CTP[ValActP[i]]
            , xylims , frecuencia )
    end
    TempInN = findall(x->x==frame,IniciosN)
    for i in TempInN
        ActivosN[ i ] = [ TN2[i][:,1] , TN2[i][:,2] ]
    end
    TempFinN = findall(x->x==frame,FinalesN)
    for i in TempFinN
        delete!(ActivosN,i)
    end
    d = 0
    for i in 1:length(ValActN)
        d+=1
        A = GraficaDinamica(c,frame,longtray,TN2[ValActN[i]]
            ,IniciosN[ValActN[i]],FinalesN[ValActN[i]],CTN[ValActN[i]]
            , xylims , frecuencia )
    end
    B = heatmap( MVcsda[:,:,frame] , colormap =:RdBu , clims = limites , tickfontsize = 8 
                , colorbar_title = "I(t)" , title = "Analisis de densidad de corriente" , 
                titlefontsize = 10 );
    if Headers != nothing
        N0s = Int( length( string( Int( FF - F1 ) ) ) )
        #println( FF, " " , F1, " " , N0s )
        NameEEG = "EEG_TOPO"
        C = Voltage_Record( frame , Headers , Canales )
        NombreFig = string( NameEEG , lpad( Int( frame - F1 + 1 ) , N0s , "0" ) , ".png" )
        DT = joinpath( ruta_pngs , NombreFig )
        Fig = load( DT )
        D = plot( Fig , axis = false )
    else
        C = nothing
        D = nothing
    end
    
    return ActivosP , ActivosN , A , B, C, D
end
#---------------------------------------------------- Crea los topogramas ----------------------------------------------------#
# Descripción: Crea las imagenes de los topogramas de las corrientes usando la función eeg_topoplot
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                                                 pesos , CoorCSDA , ruta_resultados , cachos , bines
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#         Topoplots , JLD , Plots                                   eeg_topoplot
function Crea_Topo_EEG( pesos , CoorCSDA , ruta_resultados , cachos , bines )
    NameEEG = "EEG_TOPO"    
    tamaño = floor( length( pesos ) / cachos )
    N0s = Int( length( string( Int( tamaño ) ) ) )
    NGs = length( string( cachos ) )
    Cords = CoorCSDA
    Axs = (aspect=DataAspect(),)
    Cont = (color=:black, linewidth=0.1)  
    for i in bines
        ruta_pngs = joinpath( ruta_resultados , "Temporal" )
        mkpath( ruta_pngs )
        start = Int( ( tamaño * ( i - 1) + 1 ) )
        finish = Int( ( tamaño * i ) )
        for ν in 1:tamaño
            NombreFig = string( NameEEG , lpad( Int( ν ) , N0s , "0" ) , ".png" )
            DirToSave = joinpath( ruta_pngs , NombreFig )
            x = eeg_topoplot(pesos[Int(start + ν - 1)]; positions=Cords, contours=Cont, axis=Axs)
            save( DirToSave , x ; px_per_unit = 2 )
        end
        anima = @animate for ξ in 1:tamaño
            NombreFig = string( NameEEG , lpad( Int( ξ ) , N0s , "0" ) , ".png" )
            DT = joinpath( ruta_pngs , NombreFig )
            Fig = load( DT )
            A = plot( Fig , axis = false )
        end
        ruta_gifs = joinpath( ruta_resultados * "_GIFs" );
        NomGifEEG = string( NameEEG , lpad( i , NGs , "0" ) , ".gif" )
        NomGifEEG =  joinpath( ruta_gifs , NomGifEEG )
        gif(anima, NomGifEEG , fps = 60 );
    end
    Listo = "El proceso terminó"
    return Listo
end
#--------------------------------------------------- Carga los topogramas ----------------------------------------------------#
# Descripción: Crea las imagenes de los topogramas de las corrientes usando la función eeg_topoplot
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                                                 pesos , CoorCSDA , ruta_resultados , cachos , bines
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#               JLD , Plots                             
function Carga_Topo_EEG( pesos , CoorCSDA , ruta_resultados , cachos , bines )
    # bines, tamaño, NameEEG, N0s, ruta_pngs
    for i in bines
        anima = @animate for ξ in 1:tamaño
            
            NombreFig = string( NameEEG , lpad( Int( ξ ) , N0s , "0" ) , ".png" )
            DT = joinpath( ruta_pngs , NombreFig )
            Fig = load( DT )
            A = plot( Fig , axis = false )
        end
        ruta_gifs = joinpath( ruta_resultados * "_GIFs" );
        NomGifEEG = string( NameEEG , lpad( i , NGs , "0" ) , ".gif" )
        NomGifEEG =  joinpath( ruta_gifs , NomGifEEG )
        gif(anima, NomGifEEG , fps = 60 );
    end
    Listo = "El proceso terminó"
    return Listo
end
#--------------------------------------------------- Grafica Voltaje ---------------------------------------------------------#
# Descripción: Crea las imagenes de los topogramas de las corrientes usando la función eeg_topoplot
#      Lo que nos interesa recuperar es      ///           Las variables que necesitamos
#                   Volts                                    frame , Headers , Canales
#      Las paqueterias que necesitamos       ///          Las funciones que necesitamos
#                   Plots       
function Voltage_Record( frame , Headers , Canales )
    Δ = 0
    Volts = nothing
    for Head in Headers
        largo = length( Canales[ Headers[ 1 ] ] )
        if Δ == 0
            if frame > 150 && frame < largo - 150
                Volts = plot( Canales[ Head ][ frame - 150 : frame + 150 ] , label = Head
                          , legend =:outerright )
                
            else
                if frame <= 150
                    Graficable = vcat( Float32.(zeros(151-frame)) , Canales[ Head ][ 1: frame+150] )
                    Volts = plot( Graficable , label = Head , legend =:outerright )
                else
                    Graficable = vcat( Canales[ Head ][frame-150:end],zeros( 150 -(largo - frame)) )
                    Volts = plot( Graficable , label = Head , legend =:outerright )
                end
            end
        else
            if frame > 150 && frame < largo - 150
                Volts = plot!( Canales[ Head ][ frame - 150 : frame + 150 ] .+ Δ , label = Head
                              , legend =:outerright , foreground_color_text=:white 
                              , legend_font_pointsize = 7 )
            else
                if frame <= 150
                    Graficable = vcat( Float32.(zeros(151-frame)) , Canales[ Head ][ 1: frame+150] )
                    Volts = plot!( Graficable .+ Δ, label = Head , legend =:outerright , 
                                   foreground_color_text=:white , legend_font_pointsize = 7 )
                else
                    Graficable = vcat( Canales[ Head ][frame-150:end],zeros( 150 -(largo - frame)) )
                    Volts = plot!( Graficable .+ Δ, label = Head , legend =:outerright , 
                                   foreground_color_text=:white , legend_font_pointsize = 7 )
                end
            end
        end
        Δ = Δ - 100
    end
    plot!( [ 150 , 150 ] , [ 100 , -1900 ] , color = :black , linestyle = :solid , label = frame )
    return Volts
end
end